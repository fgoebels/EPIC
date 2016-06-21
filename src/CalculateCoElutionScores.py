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
from sklearn import svm
from sklearn import datasets
from sklearn import metrics
import random
import GoldStandard as GS

subfldr = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"TCSS")))
sys.path.append(subfldr)

from main import load_semantic_similarity, calculate_semantic_similarity

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
class Apex(object):
	def __init__(self):
		self.name="apex"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))


	def setAll(self, elution):
		self.apex_array = np.argmax(elution, axis=1)
		self.shape = (len(self.apex_array),len(self.apex_array))

	def calculateScore(self,a,b):
		self.setAll(np.array([a,b]))
		return self.apex_scores_toarray_fast()[0][1]

	def __getitem__(self, index):
		 return int(self.apex_array[index[0]] == self.apex_array[index[1]])

	def apex_scores_toarray_fast(smat):
		dmaxes = defaultdict(set)
		for row, mx in enumerate(smat.apex_array):
			dmaxes[mx].add(row)
		arr = np.zeros(smat.shape)
		for mx,rows in dmaxes.items():
			for r1,r2 in itertools.permutations(rows,2):
				arr[r1,r2] = 1
		return arr

# @ author Florian Goebels
# returns bayes correlation for two proteins
# reference on how the score is calculated is here: http://www.perkinslab.ca/sites/perkinslab.ca/files/Bayes_Corr.R
class Bayes_corr:
	def __init__(self, corNum = 1):
		r=robjects.r
		r.source("src/Bayes_Corr.R")
		corFunName = "Bayes_Corr_Prior%i" % (corNum)
		self.name = corFunName
		self.corFun = robjects.r[corFunName]
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a,b):
		dims = np.matrix([a,b]).shape
		r_cro_mat_object = robjects.r.matrix(robjects.IntVector(np.append(a, b)), nrow=dims[1], ncol = dims[0])
		return self.corFun(r_cro_mat_object.transpose())[1]

# @ author Florian Goebels
# returns weighted cross correlation similarity score which is based on previously published methods Pierre 2011 et al.
class Wcc:
	def __init__(self):
		self.name="wcc"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a,b):
		rpackages.importr('wccsom')
		r_wcc = robjects.r['wcc']
		return r_wcc(robjects.FloatVector(a), robjects.FloatVector(b), 20)[0]


# @ author Florian Goebels
# returns travor correlation which is pearson correlation plus poisson noise to remove the influence of low counts
class Poisson:
	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.noiseMats = []

	def getScores(self, a, b, elutionData):
		if len(self.noiseMats)> 0 and elutionData.elutionMat.shape != self.noiseMats[0].shape:
			self.noiseMats = []
		
		if len(self.noiseMats)<self.repeat:
			for i in range(self.repeat):
				self.noiseMats.append(self.makenoisyMat(elutionData.elutionMat))
		return (elutionData.getProtIndex(a), elutionData.getProtIndex(b))

	def calculateScore(self, a,b):
		out = []
		for mat in self.noiseMats:
			profile_a = mat[a].getA1()
			profile_b = mat[b].getA1()
			mat_correlation = scipy.stats.pearsonr(profile_a, profile_b)[0]
			out.append(mat_correlation)
		return sum(out)/len(out)

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
class MutualInformation():
	def __init__(self, minCounts = 2):
		self.name="MI"
		self.minCounts = minCounts
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a, b):
		numFracs = len(a)
		(a_upper, a_lower) = self.getFracs(a, self.minCounts)
		(b_upper, b_lower) = self.getFracs(b, self.minCounts)
		entropy_a = self.bin_entropy(len(a_upper)/numFracs)
		entropy_b = self.bin_entropy(len(b_upper)/numFracs)
		joint_probs = np.array(map(lambda x: len(x)/numFracs, [a_upper&b_upper, a_upper&b_lower, a_lower&b_upper, a_lower&b_lower]))
		joint_entropy_a_b = self.entropy(joint_probs, 2)
		mutual_information =  entropy_a  + entropy_b - joint_entropy_a_b
		return mutual_information

	def bin_entropy(self, p):
		return self.entropy(np.array([p,1-p]))

	def entropy(self, probs, base=0):
		if base ==0: base = len(probs)
		tmp_probs = probs
		tmp_probs[tmp_probs==0] = 1
		return -sum(probs*map(lambda x: math.log(x,base), tmp_probs))

	def getFracs(self, a, cutoff):
		upper = set([i for i,v in enumerate(a) if v > cutoff])
		lower = set([i for i,v in enumerate(a) if v <= cutoff])
		return (upper, lower)
# @ author Florian Goebels
# calculates Jaccard overlap score which is as follows
# Jaccard(x,y) = sum(x!=1 and y!=1)/(sum(x!=1) + sum(y!=1))
class Jaccard():
        def __init__(self):
                self.name="Jaccard"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def calculateScore(self, a,b):
		j11 = 0
		j01 = 0
		for i in range(len(a)):
			if a[i] == 0 and b[i]==0: continue
			if a[i] > 0 and b[i] > 0 :
				j11 += 1
			else:
				j01 += 1
		if j11+j01 > 0:
			return j11/(j11+j01)
		else:
			return 0
# @ author Florian Goebels
# return Pearson correlation of two proteins
class Pearson:
	def __init__(self):
		self.name = "Pearson"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))


	def calculateScore(self, a,b):
		score = scipy.stats.pearsonr(a, b)[0]
		if math.isnan(score): return 0.0
		return scipy.stats.pearsonr(a, b)[0]

# @ author Florian Goebels
# returns Euclidean distance of two proteins
class Euclidiean:
	def __init__(self):
		self.name = "Euclidiean"

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a, normed=True), elutionData.getElution(b, normed=True))

	def calculateScore(self, a,b):
		return 1-distance.euclidean(a,b)

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
		numFeature_in_merge = len(toMerge.scores[toMerge.scores.keys()[0]])
		numFeature_in_self = len(self.scores[self.scores.keys()[0]])
		for edge in toMerge.scores:
			if edge in self.scores:
				self.scores[edge] = np.append(self.scores[edge], toMerge.scores[edge])
			else:
				self.scores[edge] = np.append(np.array([0]*numFeature_in_self), toMerge.scores[edge])
		for edge in self.scores:
			if edge not in toMerge.scores:
				self.scores[edge] = np.append(self.scores[edge], np.array([0]*numFeature_in_merge))


	# @author: Florian Goebels
	# returns a list with all possible pairs of protein interactions for a given ElutionData object
	def getAllPairs(self):
		allprots = self.elutionData.prot2Index.keys()
		allPPIs = set([])
		for i in range(len(allprots)):
			for j in range(i+1, len(allprots)):
				protA = allprots[i]
				protB = allprots[j]
				if protA == protB: continue
				if protA > protB:
					allPPIs.add((protA, protB, "?"))
				else:
					allPPIs.add((protB, protA, "?"))
					
		return list(allPPIs)

	# @author: Florian Goebels
	# calculates all given scores for all possible interactions of a given ElutionData object
	# @Param:
	#		scores a list of scores to calculate
	def calculateAllPairs(self, scores):
		allprots = self.getAllPairs()
		self.calculateAllScores(scores, allprots)		

	# @author: Florian Goebels
	# create co elution table for a given list of scores and ppis
	# @Param:	
	#		scoreTypes	a list of co elutoin score objects
	#		PPIs		a list of ppis for which all scores in scoreTypes schoukd be calculated
	def calculateAllScores(self, scoreTypes, PPIs):
		for scoreType in scoreTypes:
			self.header = np.append(self.header, scoreType.name)
		for protA, protB, label in PPIs:
			if (protA, protB, label) not in self.scores: self.scores[(protA, protB, label)] = []
			if not self.elutionData.hasProt(protA) or not self.elutionData.hasProt(protB):
				self.scores[(protA, protB, label)] = np.append(self.scores[(protA, protB, label)], 0)
				continue
			for scoreType in scoreTypes:
				profileA, profileB = scoreType.getScores(protA, protB, self.elutionData)
				self.scores[(protA, protB, label)] = np.append(self.scores[(protA, protB, label)], scoreType.calculateScore(profileA, profileB))


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
		for idA, idB, label in self.scores:
			if label == "positive": targets.append(1)
			if label == "negative": targets.append(0)
			if label == "?": targets.append("?")
			data.append(self.scores[(idA, idB, label)])
		return np.nan_to_num(np.array(data)), np.array(targets)

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
	print len(reference.goldstandard_positive)
	print len(reference.goldstandard_negative)
	print len(reference.goldstandard)
	return reference.goldstandard

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
# filler main function 
def main():
	refF, elutionFiles, outF = sys.argv[1:]
	elutionFH = open(elutionFiles)
	scoreCals = []
	elutionDatas = []
	
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)

	reference = createGoldStandard(elutionDatas)
	out = []
	global_all_scoreCalc = []
#	change here which features you want to use for calculating co elution scores
#	for score in [Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), MutualInformation(2)]:
#	for score in [Bayes_corr(), Bayes_corr(2), Bayes_corr(3)]:
	for score in [MutualInformation(2), Euclidiean(), Pearson(), Wcc(), Apex(), Jaccard(), Poisson(10), Bayes_corr(3)]:
		for elutionD in elutionDatas:
			scoreCalc = CalculateCoElutionScores(elutionD)
			scoreCalc.calculateAllScores([score], reference)
			scoreCals.append(scoreCalc)
			global_all_scoreCalc.append(scoreCalc)

# This code plot each feature as PR curve individually
#		all_scoreCalc = scoreCals[0]
#		for i in range(1,len(scoreCals)):
#			all_scoreCalc.mergeScoreCalc(scoreCals[i])
#		data, targets = all_scoreCalc.toSklearnData()
#		clf = CLF_Wrapper(data, targets)
#		print clf.getValScores()
#		precision, recall, _ =  clf.getPRcurve()
#		out.append((score.name, precision, recall))
#		scoreCals = []
	all_scoreCalc = global_all_scoreCalc[0]
	for i in range(1, len(global_all_scoreCalc)):
		all_scoreCalc.mergeScoreCalc(global_all_scoreCalc[i])

	data, targets = all_scoreCalc.toSklearnData()	
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
