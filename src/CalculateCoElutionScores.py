#!/usr/bin/python
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
import multiprocessing as mp
from functools import partial
import GoldStandard as GS

subfldr = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"TCSS")))
sys.path.append(subfldr)

#from main import load_semantic_similarity, calculate_semantic_similarity


# GLobal static objects used across the files
# Load R packages globally in order to avoid loading each package again in each thread

# Load script for calculating Bayes correlation globally
r=robjects.r
r.source("src/Bayes_Corr.R")
cor1 = robjects.r["Bayes_Corr_Prior1"]
cor2 = robjects.r["Bayes_Corr_Prior2"]
cor3 = robjects.r["Bayes_Corr_Prior3"]

# Packages required for calculating WCC
rpackages.importr('wccsom')
r_wcc = robjects.r['wcc']

# used in cacluating Poission score
minCounts = 2

# array for storing elution matrices with poission noise for PCC + Noise co-elution freature
noiseMats = []

# Default number of Threads is number of available cores
num_cores = mp.cpu_count()


# Helper functions for calculating co-elution scores on mutliple processors at the same time
def getScore(scoreFun, profileA, profileB, protnames):
	score = scoreFun(profileA, profileB)
	return (protnames, score)

def worker(task):
	fun , args = task
	return fun(*args)


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
		self.name = os.path.split(elutionProfileF)[-1]
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

	def clear(self):
		return True

#Fueatures that can be calcuated in parrallel have always a static helper function, since pythons multiprocessing package does not support pickeling of objects

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

	def clear(self):
		return True

	calculateScore = staticmethod(calculateScore_apex)

# @ author Florian Goebels
# returns bayes correlation for two proteins
# reference on how the score is calculated is here: http://www.perkinslab.ca/sites/perkinslab.ca/files/Bayes_Corr.R
def getRmat(a,b):
	dims = np.matrix([a, b]).shape
	r_cro_mat_object = robjects.r.matrix(robjects.IntVector(np.append(a, b)), nrow=dims[1], ncol=dims[0])
	r_cro_mat = r_cro_mat_object.transpose()
	return r_cro_mat

def calculateScoreBayes1(a,b):
	r_cro_mat = getRmat(a,b)
	global cor1
	return cor1(r_cro_mat)[1]

def calculateScoreBayes2(a,b):
	r_cro_mat = getRmat(a,b)
	global cor2
	return cor2(r_cro_mat)[1]

def calculateScoreBayes3(a,b):
	r_cro_mat = getRmat(a,b)
	global cor3
	return cor3(r_cro_mat)[1]

class Bayes1:
	def __init__(self):
		self.name = "Bayes1"
		self.parallel = True
	
	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def clear(self):
		return True

	calculateScore = staticmethod(calculateScoreBayes1)

class Bayes2:
	def __init__(self):
		self.name = "Bayes2"
		self.parallel = True

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def clear(self):
		return True

	calculateScore = staticmethod(calculateScoreBayes2)

class Bayes3:
	def __init__(self):
		self.name = "Bayes3"
		self.parallel = True

	def getScores(self, a, b, elutionData):
		return (elutionData.getElution(a), elutionData.getElution(b))

	def clear(self):
		return True

	calculateScore = staticmethod(calculateScoreBayes3)

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

	def clear(self):
		return True

# @ author Florian Goebels
# returns travor correlation which is pearson correlation plus poisson noise to remove the influence of low counts
def calculateScore_PCCPN(a,b):
	global noiseMats
	corrs = []
	for mat in noiseMats:
			corrs.append(scipy.stats.pearsonr(mat[a].getA1(), mat[b].getA1())[0])
	return sum(corrs)/len(corrs)

class Poisson:
	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.noiseMats = []
		self.parallel = True

	def getScores(self, a, b, elutionData):
		global noiseMats
		if len(noiseMats) == 0:
			self.setNoiseMats(elutionData)

		indexA = elutionData.getProtIndex(a)
		indexB = elutionData.getProtIndex(b)
		return (indexA, indexB)

	calculateScore = staticmethod(calculateScore_PCCPN)

	def setNoiseMats(self, elutionData):
		global noiseMats
		for i in range(self.repeat):
			noiseMats.append(self.makenoisyMat(elutionData.elutionMat))

	def makenoisyMat(self, mat):
		global noiseMats
		M = mat.shape[1]
		C = mat + 1/M
		poisson_mat = np.matrix(np.zeros(C.shape))
		for i in range(C.shape[0]):
			for j in range(M):
				poisson_mat[i,j] = np.random.poisson(C[i,j])
		poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
		return poisson_mat

	def clear(self):
		global noiseMats
		noiseMats = []

# @ author Florian Goebels
# returns Mutual Information of two proteins
# mutual information is based on entropy calculation MI(x,y) = H(x) + H(y) - H(x,y) 
# wehre H means entropy
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


	def clear(self):
		return True

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


	def clear(self):
		self.non_zero_fracs_for_prot = {}
		return True

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

	def clear(self):
		return True

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

	def clear(self):
		return True

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
	
	def initFromFile(self, scoreF, reference):
		scoreFH = open(scoreF)
		self.header = scoreFH.readline().rstrip()
		for line in scoreFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			idA, idB = sorted(linesplit[0:2])
			edge_scores = map(float, linesplit[2:])
			label = "?"
			for label_type in ["positive", "positive"]:
				if (idA, idB, label_type) in reference:
					label = label_type
			self.scores[(idA, idB, label)] = edge_scores

	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects unto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined tiwth self object
	def merge_singe_ScoreCalc(self, toMerge):
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
		self.header.extend(toMerge.header[2:])

	def filter_predictions_coelutionscore(self, cutoff = 0.5, cols = ""):
		if cols == "": cols = range(len(self.header)-2)
		todel = set([])
		for edge in self.scores:
			protA, protB, lab = edge
			if lab != "?": continue
			hasScore = False
			scores = self.scores[edge]
			for col in cols:
				if scores[col]>cutoff:
					hasScore = True
					break
			if not hasScore: todel.add(edge)
		for edge in todel:
				del self.scores[edge]
		print "removed a total of %i with coelution scores lower that 0.5" % (len(todel))

	def merge_list_ScoreCalc(self, toMerge):
			for sCalc in toMerge:
					self.merge_singe_ScoreCalc(sCalc)

	def mergeScoreCalc(self, toMerge):
			if isinstance(toMerge, list):
				self.merge_list_ScoreCalc(toMerge)
			else:
				self.merge_singe_ScoreCalc(toMerge)

	# @ author Florian Goebels
	# A filter for removing all possible protein pairs when predicting the network form elution data
	# here we decided to remove all candidate ppis with an Jaccard score of 0, since those interaction do not co-elute in any observed fraction
	def filter_interactions_Jaccard(self, ppis):
		print "start filtering"
		out = set([])
		thisJaccard = Jaccard()
		for protA, protB, label in ppis:
			if label != "?":
				out.add(tuple([protA, protB, label]))
			profileA, profileB = thisJaccard.getScores(protA, protB, self.elutionData)
			if not profileA.isdisjoint(profileB): out.add(tuple([protA, protB, label]))
		return out

	# @author: Florian Goebels
	# returns a list with all possible pairs of protein interactions for a given ElutionData object
	def getAllPairs(self, ref, filter = True):
		allprots = self.elutionData.prot2Index.keys()
		allPPIs = set([])
		print "Making all pairs"
		for i in range(len(allprots)):
			for j in range(i+1, len(allprots)):
				protA = allprots[i]
				protB = allprots[j]
				if protA == protB: continue
				protA, protB = sorted([protA, protB])
				label = "?"
				for label_type in set(["positive", "negative"]):
					if (protA, protB, label_type) in ref: label = label_type
				allPPIs.add((protA, protB, label))
		print "Num PPIs: %i" % (len(allPPIs))
		if filter:
			allPPIs = self.filter_interactions_Jaccard(allPPIs)
			print "After filtering Jaccard: %i num pairs" % (len(allPPIs))
		return allPPIs

	# @author: Florian Goebels
	# calculates all given scores for all possible interactions of a given ElutionData object
	# @Param:
	#		scores a list of scores to calculate
	def calculateAllPairs(self, scores, ref, filter):
		allprots = self.getAllPairs(ref, filter)
		self.calculateScores(scores, allprots)

	# @author: Florian Goebels
	# create co elution table for a given list of scores and ppis
	# @Param:	
	#		scoreTypes	a list of co elutoin score objects
	#		PPIs		a list of ppis for which all scores in scoreTypes schoukd be calculated
	def calculateScores(self, scoreTypes, PPIs):
		for st in scoreTypes:
			thisst = copy.deepcopy(st) # copy score object beacuse some objects store information internally and by copying them avoid having data stored when calculating the score for different data sets in the same run
			name = self.elutionData.name + "." +  st.name
			print "Calculating %s" % (name)
			global noiseMats
			print "Number of noise mpats for pcc+n:%i" % (len(noiseMats))
			numScores = 1
			if isinstance(name, list):
				self.header.extend(name)
				numScores = len(name)
			else:
				self.header.append(name)
			tasks = []
			for protA, protB, label in PPIs:
				if (protA, protB, label) not in self.scores: self.scores[(protA, protB, label)] = []
				if not self.elutionData.hasProt(protA) or not self.elutionData.hasProt(protB):
					self.scores[(protA, protB, label)].extend([0]*numScores)
					continue
				profileA, profileB = thisst.getScores(protA, protB, self.elutionData)
				task = (getScore, (thisst.calculateScore, profileA, profileB, (protA, protB, label)))
				tasks.append(task)
			if thisst.parallel == True:
				global num_cores 
				p = mp.Pool(num_cores)
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
			print "Number of noise mpats for pcc+n:%i" % (len(noiseMats))
			st.clear()
			print "Number of noise mpats for pcc+n:%i" % (len(noiseMats))

	def calculate_coelutionDatas(self, elutionDatas, scores, reference, predict, filter):
		for elutionData in elutionDatas:
			print "Processing: %s" % (elutionData.name)
			scoreCalc = CalculateCoElutionScores(elutionData)
			if predict:
				scoreCalc.calculateAllPairs(scores, reference, filter)
			else:
				scoreCalc.calculateScores(scores, reference)
			self.mergeScoreCalc(scoreCalc)

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
	def toSklearnData(self, labels = ["positive", "negative", "?"]):
		data = []
		targets = []
		ids = []
		for idA, idB, label in self.scores:
			if label not in labels: continue
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
	def __init__(self, data, targets, forest=False, folds =10):
		self.data = data
		self.folds = folds
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		self.eval_preds = ""
		self.eval_targets = ""
		if forest:
			print ("using Random forest")
			self.clf = RandomForestClassifier(n_estimators=400, n_jobs=num_cores)
		else:	
			print ("Using SVM")
			self.clf =  svm.SVC(kernel="linear", probability=True)
		self.clf.fit(self.data, self.targets)

	# @author Florian Goebels
	def geteval(self):
		all_probas = []
		all_targets = []
		for train, test in StratifiedKFold(self.targets, self.folds):
			probas = self.clf.fit(self.data[train], self.targets[train]).predict_proba(self.data[test])
			all_probas.extend(probas[:,1]) # make sure that 1 is positive class in binarizied class vector
			all_targets.extend(self.targets[test])
		self.eval_preds  = all_probas
		self.eval_targets = all_targets

	# @author: Florian Goebels
	# calculates precision, recall, fmeasure, auc_pr, auc_roc for n fold cross validation
	# @Param:
	#		folds number of folds default 10
	def getValScores(self):
		global num_cores
		folds = StratifiedKFold(self.targets, self.folds)
		preds = cross_val_predict(self.clf, self.data, self.targets, cv=folds, n_jobs=num_cores)
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
		if self.eval_preds == "":
			self.geteval()
		return precision_recall_curve(self.eval_targets, self.eval_preds)

	def getROCcurve(self, folds=10):
		if self.eval_preds == "":
			self.geteval()
		return roc_curve(self.eval_targets, self.eval_preds)

	# @author: Florian Goebels
	# @Param:
	#		toPred matric where each row is a data point and predicts interaction propability for a given set
	#		note trainer needs to be already trained to be ablte to predict
	def predict(self, toPred):
		probas = self.clf.predict_proba(toPred)
		return probas[:,1]

# @author: Florian Goebels
# makes precision recall plot for mutliple rp ccurves
# @Param:
#	prcurves list of tuples with (name, precision, recall) which should be plotted
#	outF Pdf file location for the created plot
def plotCurves(curves, outF, xlab, ylab):
	plt.clf()
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	cols = ['b', 'r', 'c', 'm', 'y', 'k', '#0000ff', '#005D00', '#00FF00', '#8B5D00', '#FF6600', '#FF0066', '#5c5c8a']
	for (name, y, x) in curves:
		if name != "":
			plt.plot(x, y, label=name, color = cols.pop())
		else:
			plt.plot(x, y, color=cols.pop())
	art = []
	lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1),  ncol = 5, fontsize=8)
	art.append(lgd)
	plt.savefig(outF, additional_artists=art, bbox_inches="tight")

# @author Florian Goebels
def predictInteractions(reference, elutionDatas, scores, outF, useForest):
	scoreCalc = CalculateCoElutionScores()
	scoreCalc.calculate_coelutionDatas(elutionDatas, scores, reference, predict=True, filter=True)
	scoreCalc.filter_predictions_coelutionscore()
	print len(reference)
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(labels=set(["positive", "negative"]))
	print data_train.shape
	ids_pred, data_pred, targets_pred = scoreCalc.toSklearnData()
	print data_pred.shape
	clf = CLF_Wrapper(data_train, targets_train, forest=useForest)
	preds = clf.predict(data_pred)
	outFH = open(outF, "w")
	cutoff = 0.1
	if useForest: cutoff = 0.6
	for i in range(len(ids_pred)):
		if preds[i] >= cutoff:
			outFH.write("%s\t%f\n" % ("\t".join(ids_pred[i]), preds[i]))
	outFH.close()


# @author Florian Goebels
def main():
	refF, useForest, elutionFiles, outDir = sys.argv[1:]
	useForest = useForest == "True"

	scores = [MutualInformation(2), Bayes1(), Bayes2(), Bayes3(), Euclidiean(), Pearson(), Wcc(), Jaccard(), Poisson(10)]
#	scores = [Pearson(), Poisson(10)]

	elutionFH = open(elutionFiles)
	elutionDatas = []
	elutionProts = set([])
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())
	
	reference = GS.Goldstandard_from_reference_File(refF, found_prots=elutionProts)
#	bench_scores(reference.goldstandard, elutionDatas, scores, outDir, useForest)
#	print len(reference.goldstandard_positive)
#	print len(reference.goldstandard)
	predictInteractions(reference.goldstandard, elutionDatas, scores, outDir + ".pred.txt", useForest)

def bench_scores(reference, elutionDatas, scores, outDir, useForest=False):
	scoreCals = []
	pr_curves = []
	roc_curves = []
	cutoff_curves = []
	global_all_scoreCalc = []
	tableFH = open(outDir + ".eval.txt", "w")
	tableFH.write("Name\tPrecision\tRecall\tF-measure\tau_PR\tau_ROC\n")
	for score in scores:
		print(score.name)
		scoreCalc = CalculateCoElutionScores()
		scoreCalc.calculate_coelutionDatas(elutionDatas, [score], reference, predict=False, filter=False)
		_, data, targets = scoreCalc.toSklearnData()
		print("starting ml")
		clf = CLF_Wrapper(data, targets, forest=useForest, folds=10)
		tableFH.write("%s\t%s\n" % (score.name, "\t".join(map(str, clf.getValScores()))))
		precision, recall, _ =  clf.getPRcurve()
		pr_curves.append((score.name, precision, recall))

		fpr, tpr, _ = clf.getROCcurve()
		roc_curves.append((score.name, tpr, fpr))
		print ("done with ml")
		scoreCals.append(scoreCalc)

	print ("starting combined ml")
	all_scoreCalc = CalculateCoElutionScores()
	all_scoreCalc.merge_list_ScoreCalc(scoreCals)
	_, data, targets = all_scoreCalc.toSklearnData()	
	clf = CLF_Wrapper(data, targets, forest=useForest, folds=10)
	tableFH.write("Combined\t%s\n" % ("\t".join(map(str, clf.getValScores()))))
	tableFH.close()
	precision, recall, threshold =  clf.getPRcurve()
	pr_curves.append(("combined", precision, recall))
	plotCurves(pr_curves, outDir + ".pr.pdf", "Recall", "Precision")

	dataFH = open(outDir + ".scores.txt", "w")
	dataFH.write(all_scoreCalc.toTable())
	dataFH.close()

	threshold = np.append(threshold, 1)
	fpr, tpr, _ = clf.getROCcurve()
	roc_curves.append(("combined", tpr, fpr))
	plotCurves(roc_curves, outDir + ".roc.pdf", "False Positive rate", "True Positive Rate")
	
	cutoff_curves.append(("Precision", precision, threshold))
	cutoff_curves.append(("Recall", recall, threshold))
	plotCurves(cutoff_curves, outDir + ".cutoff.pdf", "Cutoff", "Evaluation metric score")
		


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
