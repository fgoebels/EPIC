#!/usr/bin/env python
from __future__ import division
import numpy as np
import scipy.stats
import utils 
import sys
import math
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from scipy.spatial import distance
import itertools
import operator
import copy
from sklearn import cross_validation
from sklearn.cross_validation import StratifiedKFold
from collections import defaultdict
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.cross_validation import train_test_split, cross_val_predict
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import inspect
import os
import time
from sklearn import svm
from sklearn import datasets
from sklearn import metrics
import random
#import pathos.multiprocessing as mp
#import pathos.pp as pp
import multiprocessing as mp
from functools import partial
import GoldStandard as GS

subfldr = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"TCSS")))
sys.path.append(subfldr)

from main import load_semantic_similarity, calculate_semantic_similarity

r=robjects.r
r.source("src/Bayes_Corr.R")
cor1 = robjects.r["Bayes_Corr_Prior1"]
cor2 = robjects.r["Bayes_Corr_Prior2"]
cor3 = robjects.r["Bayes_Corr_Prior3"]

rpackages.importr('wccsom')
r_wcc = robjects.r['wcc']



#Contacts:
#	Florian Goebels - florian.goebels@gmail.com
#
# @author: Florian Goebels
# This class stores and maintains co elution data
# It stores both the raw and normalizied co-elution data
class ElutionData():
	# @author: Florian Goebels
	# Class constructor
	# @Param:
	#		elutionProfileF flat file containing co elution data where each row represents a protein and collums contain spectral counts
	# @Class objects: 
	#		self.elutionMat numpy matric cnsisting of only the counts 
	#		self.normedElutionMat numpy matric containing normalized spectral counts 
	#		self.prot2Index dictonary mapping protein to it's respective row index in self.elutionMat and self.normedElutionMat
	def __init__(self, elutionProfileF):
		self.elutionMat, self.prot2Index  = self.loadElutionData(elutionProfileF)
		self.normedElutionMat = normalize_fracs(self.elutionMat)
		self.elutionMat = np.array(self.elutionMat)
		

	# @author: Florian Goebels
	# this methods load elution data stored as a flat file and returns the read in matrix and an index pointing each protein to it's rwo index
	# @Param:
	#		elutionProfileF elution data as flat file, tab separated
	def loadElutionData(self, elutionProfileF):
		elutionProfileFH = open(elutionProfileF)
		elutionProfileFH.readline()
		elutionProfile = {}
		i = 0
		elutionMat = []
		prot2Index = {}
		for line in elutionProfileFH:
			line = line.rstrip()
			line = line.split("\t")
			protID = line[0]
			counts = np.array(map(float, line[1:]))
			elutionMat.append(counts)
			prot2Index[protID] = i
			i += 1
		elutionProfileFH.close()
		return elutionMat, prot2Index

		

	# @author: Florian Goebels
	# return nurmalized co elution matrix, Normalization is based on previously published methods Pierre 2011 et al.
	def normalizeCoEulitionMat(self):
		self.elutionMat = normalize_fracs(self.elutionMat)
		

	# @author: Florian Goebels
	# returns elution profile of a given protein. If the protein has no profile the method returns NaN
	# @Param:
	#		normed if True return normalized counts, else raw counts
	def getElution(self, prot, normed=False):
		if prot not in self.prot2Index:
			return float('NaN')
		if normed: 
			return self.normedElutionMat[self.prot2Index[prot]]
		else:
			return self.elutionMat[self.prot2Index[prot]]
	

	# @author: Florian Goebels
	# returns true if there is elution profile for a protein
	# @Param:
	#		prot  the protein in question
	def hasProt(self, prot):
		return prot in self.prot2Index

	# @author: Florian Goebels
	# reutrns a protein's row index in the co-elution table
	# @Param:
	#		prot  the protein in question
	def getProtIndex(self, prot):
		return self.prot2Index[prot]

	# @author: Florian Goebels
	# helper function for printing co-elution mat. Returns string
	# @Param:
	#		co elution mat as np matrix, which to be printed
	def printMat(self, mat):
		out = "ProtiID\tFraction_" + "\tFracion_".join(map(str, range(1,240)))
		for prot in self.prot2Index:
			index = self.prot2Index[prot]
			out += "\n%s\t%s" %  (prot, "\t".join(map(str,self.normedElutionMat[index])))
		return out

# @author: Florian Goebels
# helper function for normalization
def arr_norm(arr, axis=0):
    """
    axis=0: normalize each column; 1: normalize each row
    """
    mat = np.asmatrix(arr)
    return np.asarray(np.nan_to_num(mat / np.sum(mat, axis)))

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
    if norm_cols:
        # Normalize columns first--seems correct for overall elution profile if
        # you're calculating correlation-type scores
        arr = arr_norm(arr, 0)
    if norm_rows:
        arr = arr_norm(arr, 1)
    return arr


# The following classes describe features for calculating Co elution scores
# each class needs to have two central functions: getScores and calculateScore
# whereas the first one allows the method to prepare the elution data how ever it needs and return int as two objects
# and the later one return the score for two given score objects


# @author: Florian Goebels
# Calculates functional similarity based on GO terms
class GOSim(object):

	objs = "" 
	def __init__(self, gene_anno_F, onto_F):
		self.name = ["Sim_CC", "Sim_BP", "Sim_MF"]
		self.parallel = False
#		if GOSim.objs == "":
#			GOSim.objs = load_semantic_similarity(onto_F, gene_anno_F, "C:2.4,P:3.5,F:3.3", "IEA")

	def getScores(self, a, b, elutionData):
		return (a,b)

	def calculateScore(self, a, b):
		out = []
		domain_def = {'C':'Cellular Component', 'P':'Biological Process', 'F':'Molecular Function'}
		for domain in domain_def:
			score = GOSim.objs[domain]._semantic_similarity(a, b)[0]
			if score is None: score = 0
			out.append(score)
		return out

# @ author Florian Goebels
# returns apex similarity score which is based on previously published methods Pierre 2011 et al.
# returns 0 or 1 depending if the highest count peak of two proteins are in the same fraction or not
def calculateScore_apex(a,b):
	if a == b:
		return 1
	else:
		return 0

class Apex(object):
	def __init__(self):
		self.name="apex"
		self.parallel = True

	def getScores(self, a, b, elutionData):
		return (np.argmax(elutionData.getElution(a)), np.argmax(elutionData.getElution(b)))

	calculateScore = staticmethod(calculateScore_apex)

# @ author Florian Goebels
# returns bayes correlation for two proteins
# reference on how the score is calculated is here: http://www.perkinslab.ca/sites/perkinslab.ca/files/Bayes_Corr.R

def calculateScoreBayes(a,b):
	global cor1
	global cor2
	global cor3
	dims = np.matrix([a,b]).shape
	r_cro_mat_object = robjects.r.matrix(robjects.IntVector(np.append(a, b)), nrow=dims[1], ncol = dims[0])
	r_cro_mat = r_cro_mat_object.transpose()
	return [cor1(r_cro_mat)[1], cor2(r_cro_mat)[1], cor3(r_cro_mat)[1]]
	
class Bayes_corr:
	def __init__(self, corNum = 1):
		self.name = ["Bayes1", "Bayes2", "Bayes3"]
		self.parallel = True
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	calculateScore = staticmethod(calculateScoreBayes)

# @ author Florian Goebels
# returns weighted cross correlation similarity score which is based on previously published methods Pierre 2011 et al.
def calculateScore_wcc(a,b):
	global r_wcc
	return r_wcc(robjects.FloatVector(a), robjects.FloatVector(b), 20)[0]

class Wcc:
	def __init__(self):
		self.name="wcc"
		self.parallel = True

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	calculateScore = staticmethod(calculateScore_wcc)


# @ author Florian Goebels
# returns travor correlation which is pearson correlation plus poisson noise to remove the influence of low counts
def calculateScore_PCCPN(a,b):
	corrs = map(lambda x : scipy.stats.pearsonr(a[x], b[x])[0], range(len(a)))
	return sum(corrs)/len(corrs)

class Poisson:
	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.noiseMats = []
		self.parallel = True

	def getScores(self, a, b, elutionData):
		if len(self.noiseMats)> 0 and elutionData.elutionMat.shape != self.noiseMats[0].shape:
			self.noiseMats = []
		
		if len(self.noiseMats)<self.repeat:
			for i in range(self.repeat):
				self.noiseMats.append(self.makenoisyMat(elutionData.elutionMat))
		indexA = elutionData.getProtIndex(a)
		indexB = elutionData.getProtIndex(b)
		outA = []
		outB = []
		for mat in self.noiseMats:
			outA.append(mat[indexA].getA1())
			outB.append(mat[indexB].getA1())
		return (outA, outB)

	calculateScore = staticmethod(calculateScore_PCCPN)

	def makenoisyMat(self, mat):
		M = mat.shape[1]
		C = mat + 1/M
		poisson_mat = np.matrix(np.zeros(C.shape))
		for i in range(C.shape[0]):
			for j in range(M):
				poisson_mat[i,j] = np.random.poisson(C[i,j])
		poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
		return poisson_mat

# @ author Florian Goebels
# returns Mutual Information of two proteins
# mutual information is based on entropy calculation MI(x,y) = H(x) + H(y) - H(x,y) 
# wehre H means entropy


minCounts = 2
def calculateScore_MI( a, b):
	numFracs = len(a)
	(a_upper, a_lower) = getFracs(a, minCounts)
	(b_upper, b_lower) = getFracs(b, minCounts)
	entropy_a = bin_entropy(len(a_upper)/numFracs)
	entropy_b = bin_entropy(len(b_upper)/numFracs)
	joint_probs = np.array(map(lambda x: len(x)/numFracs, [a_upper&b_upper, a_upper&b_lower, a_lower&b_upper, a_lower&b_lower]))
	joint_entropy_a_b = entropy(joint_probs, 2)
	mutual_information =  entropy_a  + entropy_b - joint_entropy_a_b
	return mutual_information

def bin_entropy(p):
	return entropy(np.array([p,1-p]))

def entropy(probs, base=0):
	if base ==0: base = len(probs)
	tmp_probs = probs
	tmp_probs[tmp_probs==0] = 1
	return -sum(probs*map(lambda x: math.log(x,base), tmp_probs))

def getFracs( a, cutoff):
	upper = set([i for i,v in enumerate(a) if v > cutoff])
	lower = set([i for i,v in enumerate(a) if v <= cutoff])
	return (upper, lower)

class MutualInformation():
	def __init__(self, minCounts_new = 2):
		self.name="MI"
		global minCounts
		minCounts = minCounts_new
		self.parallel = True
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	calculateScore = staticmethod(calculateScore_MI)

# @ author Florian Goebels
# calculates Jaccard overlap score which is as follows
# Jaccard(x,y) = sum(x!=1 and y!=1)/(sum(x!=1) + sum(y!=1))
def calculateScore_Jaccard(a_non_zero_fracs,b_non_zero_fracs):
	a_and_b_non_zero_fracs = len(a_non_zero_fracs & b_non_zero_fracs)
	a_or_b_non_zero_fracs = len(a_non_zero_fracs | b_non_zero_fracs)
	if a_or_b_non_zero_fracs == 0:
		return 0
	else:
		return a_and_b_non_zero_fracs/a_or_b_non_zero_fracs

class Jaccard():
        def __init__(self):
                self.name="Jaccard"
		self.non_zero_fracs_for_prot = {}
		self.parallel = True

	def getScores(self, a, b, elutionData):
		for prot in [a,b]:
			if prot not in self.non_zero_fracs_for_prot:
				self.non_zero_fracs_for_prot[prot] = set(np.nonzero(elutionData.getElution(prot))[0])
		return (self.non_zero_fracs_for_prot[a],self.non_zero_fracs_for_prot[b])

	calculateScore = staticmethod(calculateScore_Jaccard)


# @ author Florian Goebels
# return Pearson correlation of two proteins


def calculateScore_Pearson(a,b):
	score = scipy.stats.pearsonr(a, b)[0]
	if math.isnan(score): return 0.0
	return scipy.stats.pearsonr(a, b)[0]


class Pearson:
	def __init__(self):
		self.name = "Pearson"
		self.parallel = True

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	
	calculateScore = staticmethod(calculateScore_Pearson)

# @ author Florian Goebels
# returns Euclidean distance of two proteins
def calculateScore_euclidean(a,b):
	return 1-distance.euclidean(a,b)

class Euclidiean:
	def __init__(self):
		self.name = "Euclidiean"
		self.parallel = True


	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, normed=True), elutionData.getElution(b, normed=True))

	calculateScore = staticmethod(calculateScore_euclidean)

def getScore(scoreFun, profileA, profileB, protnames):
	score = scoreFun(profileA, profileB)
	return (protnames, score)

def worker(task):
	fun , args = task
	return fun(*args)

# @ author Florian Goebels
# This is a helper class for calculating co elution scores for a given ElutionData object
class CalculateCoElutionScores():
	# @author: Florian Goebels
	# this method inits the object
	# @Param:
	#		elutionData (optional) specify which data needs to be processed
	def __init__(self, elutionData=""):
		self.elutionData = elutionData
		self.scores = {}
		self.header = ["ProtA","ProtB"]
		

	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects unto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined tiwth self object
	def mergeScoreCalc(self, toMerge):
		numFeature_in_merge = len(toMerge.header)-2
		numFeature_in_self = len(self.header)-2
		for edge in toMerge.scores:
			if edge in self.scores:
				self.scores[edge] = np.append(self.scores[edge], toMerge.scores[edge])
			else:
				self.scores[edge] = np.append(np.array([0]*numFeature_in_self), toMerge.scores[edge])
		for edge in self.scores:
			if edge not in toMerge.scores:
				self.scores[edge] = np.append(self.scores[edge], np.array([0]*numFeature_in_merge))

	#
	#
	@staticmethod
	def filter_interactions_Jaccard(elutionData, ppis):
		thisCS = CalculateCoElutionScores(elutionData)
		thisCS.calculateAllScores([Jaccard()], ppis)
		toDel = set([])
		for edge in thisCS.scores:
			if edge[2] != "?": continue #only remove prediction label not positive or negative label name
			jaccardScore = thisCS.scores[edge][-1]
			if jaccardScore == 0.0:
				toDel.add(edge)
		out = set(ppis) - toDel
		return list(out)	

	# @author: Florian Goebels
	# returns a list with all possible pairs of protein interactions for a given ElutionData object
	def getAllPairs(self, ref):
		allprots = self.elutionData.prot2Index.keys()
		allPPIs = set([])
		for i in range(len(allprots)):
			for j in range(i+1, len(allprots)):
				protA = allprots[i]
				protB = allprots[j]
				if protA == protB: continue
				protA, protB = sorted([protA, protB])
				label = "?"
				for label_type in ["positive", "positive"]:
					if (protA, protB, label_type) in ref: label = label_type
				allPPIs.add((protA, protB, label))
		return list(allPPIs)

	# @author: Florian Goebels
	# calculates all given scores for all possible interactions of a given ElutionData object
	# @Param:
	#		scores a list of scores to calculate
	def calculateAllPairs(self, scores, ref):
		allprots = self.getAllPairs(ref)
		allprots = CalculateCoElutionScores.filter_interactions_Jaccard(self.elutionData, allprots)
		self.calculateAllScores(scores, allprots)		

	# @author: Florian Goebels
	# create co elution table for a given list of scores and ppis
	# @Param:	
	#		scoreTypes	a list of co elutoin score objects
	#		PPIs		a list of ppis for which all scores in scoreTypes schoukd be calculated
	def calculateAllScores(self, scoreTypes, PPIs):
		for st in scoreTypes:
			name = st.name
			numScores = 1
			if isinstance(name, list):
				self.header.extend(name)
				numScores = len(name)
			else:
				self.header.append(name)
			scores = []
			tasks = []
			results = []
			for protA, protB, label in PPIs:
				if (protA, protB, label) not in self.scores: self.scores[(protA, protB, label)] = []
				if not self.elutionData.hasProt(protA) or not self.elutionData.hasProt(protB):
					
					self.scores[(protA, protB, label)].extend([0]*numScores)
					continue
				profileA, profileB = st.getScores(protA, protB, self.elutionData)
				task = (getScore, (st.calculateScore, profileA, profileB, (protA, protB, label)))
				tasks.append(task)
#			results = p.map_async(worker, tasks)
#			results.wait()
#			results = results.get()
			if st.parallel == True:
				p = mp.Pool(4)
				results = p.map(worker, tasks)
				p.close()
				p.join()
			else:
				results = map(worker, tasks)
			for ppi, score in results:
				if isinstance(score, list):
					self.scores[ppi].extend(score)
				else:
					self.scores[ppi].append(score)	

	# @author: Florian Goebels
	# prints table
	# @Param:
	#		labels print class lable or not
	def toTable(self, labels=True):
		out = "\t".join(self.header)
		if labels: out += "\tClass"
		for protA, protB, label in self.scores:
			out += "\n%s\t%s\t%s" % (protA, protB, "\t".join(map(str,self.scores[(protA, protB, label)])))
			if labels: out += "\t" + label
		return out
	

	# @author: Florian Goebels
	# prints rff table
	def toArffData(self):
		out = ["@RELATION COFrac"]
		for colname in self.header[2:]:
			out.append("@Attribute %s NUMERIC" % (colname))
		out.append("@ATTRIBUTE class {positive, negative}")
		out.append("@Data")
		for idA, idB, label in self.scores:
			tmp = ",".join(map(str,self.scores[(idA, idB, label)]))
			tmp += ",%s" % (label)
			out.append(tmp)
		return "\n".join(out)

	# @author: Florian Goebels
	# return stored co elution scores for a given data set in form that is required for sklearn to learn and predict interaction
	def toSklearnData(self):
		data = []
		targets = []
		ids = []
		for idA, idB, label in self.scores:
			if label == "positive": targets.append(1)
			if label == "negative": targets.append(0)
			if label == "?": targets.append("?")
			data.append(self.scores[(idA, idB, label)])
			ids.append(tuple([idA, idB]))
		return ids, np.nan_to_num(np.array(data)), np.array(targets)

# @ author Florian Goebels
# wrapper for machine learning
class CLF_Wrapper:
	# @author: Florian Goebels
	# class initializer, supports both random forest and svm
	# @Param:
	#		data matrix with features where each row is a data point
	#		targets list with class lables
	# 		forest if true ml is random forst, if false ml is svm
	def __init__(self, data, targets, forest=False):
		self.data = data
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		if forest:
			self.clf = RandomForestClassifier(n_estimators=100)
		else:	
			self.clf = svm.SVC(kernel="linear", probability=True)
		self.clf.fit(data, targets)
		

	# @author: Florian Goebels
	# do 10 fold cross validation returns predicted class lables
	# @Param:
	#		folds number of folds default 10
	def kFoldCV(self, folds=10):
		folds = StratifiedKFold(self.targets, folds)
		return cross_val_predict(self.clf, self.data, self.targets, cv=folds)


	# @author: Florian Goebels
	# calculates precision, recall, fmeasure, auc_pr, auc_roc for n fold cross validation
	# @Param:
	#		folds number of folds default 10
	def getValScores(self, folds=10):
		preds = self.kFoldCV(folds)
		precision = metrics.precision_score(self.targets, preds, average=None)[1]
		recall = metrics.recall_score(self.targets, preds, average=None)[1]
		fmeasure = metrics.f1_score(self.targets, preds, average=None)[1]
		auc_pr = average_precision_score(self.targets, preds)
		auc_roc = roc_auc_score(self.targets, preds) 
		return [precision, recall, fmeasure, auc_pr, auc_roc]


	# @author: Florian Goebels
	# calculated PR curve using 10 fold cross validation
	# returns three arrays: precision, recall, thresholds
	# @Param:
	#		folds number of folds default 10
	def getPRcurve(self, folds=10):
		all_probas = []
		all_targets = []
		for train, test in StratifiedKFold(self.targets, folds):
			probas = self.clf.fit(self.data[train], self.targets[train]).predict_proba(self.data[test])
			all_probas.extend(probas[:,1]) # make sure that 1 is positive class in binarizied class vector
			all_targets.extend(self.targets[test])
		return precision_recall_curve(all_targets, all_probas)
#		return roc_curve(all_targets, all_probas)
	

	# @author: Florian Goebels
	# @Param:
	#		toPred matric where each row is a data point and predicts interaction propability for a given set
	#		note trainer needs to be already trained to be ablte to predict
	def predict(self, toPred):
		return self.clf.predict_proba(toPred)

# @author: Florian Goebels
# Create gold standard data set from CORUM using the GoldStandard helper class
# @Param:
#		ratio the size realtion ship between positive and negative training data points: |negative| = ratio * |positive|
#		target pecies to predict in given in taxid here we use worm i.e. taxid 6239 (mouse 10090, Human 9606) 
def createGoldStandard(elutionData, ratio=5, targetSpecies = "6239"):
	protsWithElutionProfile = set([])
	for data in elutionData:
		protsWithElutionProfile = protsWithElutionProfile | set(data.prot2Index.keys())

	inparanoid = GS.Inparanoid(targetSpecies, foundProts = protsWithElutionProfile)
	gos = GS.QuickGO(targetSpecies)
	corum = GS.CORUM()
	reference = GS.Goldstandard_from_CORUM(corum, inparanoid, ratio=5, found_prots = protsWithElutionProfile)
	return reference.goldstandard

def createGoldStandard_alt(refF, elutionData, ratio=5):
        positive = set([])
        negative = set([])
        protsWithElutionProfile = set([])
        for data in elutionData:
                protsWithElutionProfile = protsWithElutionProfile | set(data.prot2Index.keys())
        reference = set([])
        goldstandardFH = open(refF)
        for line in goldstandardFH:
                line = line.rstrip()
                (tmpidA, tmpidB, label) = line.split("\t")
                if tmpidA not in protsWithElutionProfile or tmpidB not in protsWithElutionProfile: continue
                idA, idB = sorted([tmpidA, tmpidB])
                if label == "positive": positive.add((idA, idB, label))
                if label == "negative": negative.add((idA, idB, label))

        if len(positive)*ratio>len(negative):
                print "Warning: not enough negative data points in reference to create desired ratio"
                return positive | negative

        print len(positive)
        print len(negative)
        reference = positive | set(random.sample(negative, len(positive)*ratio))
        return reference

# @author: Florian Goebels
# makes precision recall plot for mutliple rp ccurves
# @Param:
#	prcurves list of tuples with (name, precision, recall) which should be plotted
#	outF Pdf file location for the created plot
def plotPRcurve(prcurves, outF):
	plt.clf()
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	cols = ['b', 'r', 'c', 'm', 'y', 'k', '#0000ff', '#005D00', '#00FF00', '#8B5D00', '#FF6600', '#FF0066', '#5c5c8a']
	for (name, precision, recall) in prcurves:
		plt.plot(recall, precision, label=name, color = cols.pop())	
	plt.legend(loc="upper right", ncol = 5, fontsize=8)
	plt.savefig(outF)

# @author Florian Goebels
# 
def predictInteractions(taxid, elutionFile, outF):
	elutionData = ElutionData(elutionFile)
	reference = createGoldStandard([elutionData])
	scoreCalc = CalculateCoElutionScores(elutionData)
	scoreCalc.calculateAllPairs([MutualInformation(2), Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), Poisson(10), Bayes_corr()], reference)
#	scoreCalc.calculateAllPairs([Apex()], reference)
	ids, data, targets = scoreCalc.toSklearnData()
#	clf = CLF_Wrapper(data, targets)
#	preds = clf.predict(data)
	outFH = open(outF, "w")
	outFH.write(scoreCalc.toTable())
	outFH.close()
	
# @author Florian Goebels
# filler main function 
def main():
#	taxid, elutionF, outF = sys.argv[1:]
#	predictInteractions(taxid, elutionF, outF)
	refF, elutionFiles, outF = sys.argv[1:]
	bench_scores(refF, elutionFiles, outF)

def bench_scores(refF, elutionFiles, outF):
	elutionFH = open(elutionFiles)
	scoreCals = []
	elutionDatas = []
	
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)

	reference = createGoldStandard_alt(refF, elutionDatas)
	out = []
	global_all_scoreCalc = []
#	change here which features you want to use for calculating co elution scores
#	for score in [Poisson(10)]:
	for score in [MutualInformation(2), Bayes_corr(), Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), Poisson(10)]:
		print score.name
		for elutionD in elutionDatas:
			scoreCalc = CalculateCoElutionScores(elutionD)
			scoreCalc.calculateAllScores([score], reference)
			scoreCals.append(scoreCalc)
			global_all_scoreCalc.append(scoreCalc)
		
# This code plot each feature as PR curve individually
		all_scoreCalc = scoreCals[0]
		for i in range(1,len(scoreCals)):
			all_scoreCalc.mergeScoreCalc(scoreCals[i])
		_, data, targets = all_scoreCalc.toSklearnData()
		print "starting ml"
		clf = CLF_Wrapper(data, targets)
		print clf.getValScores()
		precision, recall, _ =  clf.getPRcurve()
		out.append((score.name, precision, recall))
		scoreCals = []
		print "done with ml"
	all_scoreCalc = global_all_scoreCalc[0]
	print "strting combined ml"
	for i in range(1, len(global_all_scoreCalc)):
		all_scoreCalc.mergeScoreCalc(global_all_scoreCalc[i])

	_, data, targets = all_scoreCalc.toSklearnData()	
	clf = CLF_Wrapper(data, targets)
	print clf.getValScores()
	precision, recall, _ =  clf.getPRcurve()
	out.append(("combined", precision, recall))

	plotPRcurve(out, outF)
	

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pas
