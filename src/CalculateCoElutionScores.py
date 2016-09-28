#!/usr/bin/python
from __future__ import division
import numpy as np
import scipy.stats
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
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFECV
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
from ctypes import c_wchar, c_float
import logging

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
current_bayes_cor = ""

# Packages required for calculating WCC
rpackages.importr('wccsom')
r_wcc = robjects.r['wcc']

# array for storing elution matrices with poission noise for PCC + Noise co-elution freature
Poisson_cor_Mat = []

# Default number of Threads is number of available cores
num_cores = mp.cpu_count()

#Global variables for calculating scores across mutliple threads
prot2profile = {}
ppiList = []


# Helper functions for calculating co-elution scores on mutliple processors at the same time
def getScore(inputQ, outputList):
	for args in iter(inputQ.get, 'STOP'):
		global prot2profile
		ppi, score_index, fname, scoreFun = args
		protA, protB = ppi.split("\t")
		if protA in prot2profile[fname] and protB in prot2profile[fname]:
			profileA = prot2profile[fname][protA]
			profileB = prot2profile[fname][protB]
			score = scoreFun(profileA, profileB)
			outputList.put((score_index, score))
		else:
			outputList.put((score_index, 0.0))
		inputQ.task_done()

def worker(input, output, results):
	for func, args in iter(input.get, 'STOP'):
		ppi, score = func(*args)
		if float(score) > 0.5: results.append(ppi)
		output.put("%s\t%s" % (ppi, score))
		input.task_done()

def listener(output, fh):
	for line in iter(output.get, 'STOP'):
		print >> fh, line
		fh.flush()
		output.task_done()

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
		i = 0
		elutionMat = []
		prot2Index = {}
		for line in elutionProfileFH:
			line = line.rstrip()
			line = line.split("\t")
			protID = line[0]
			counts = map(float, line[1:])
			elutionMat.append(counts)
			prot2Index[protID] = i
			i += 1
		elutionProfileFH.close()
		elutionMat = np.nan_to_num(np.matrix(elutionMat))
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
			return None
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
		if prot in self.prot2Index:
			return self.prot2Index[prot]
		else:
			return None

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

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = np.argmax(elutionData.getElution(prot))

	def clear(self):
		global prot2profile
		prot2profile = {}
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

def calculateScoreBayes(a,b):
	r_cro_mat = getRmat(a,b)
	global current_bayes_cor
	return  current_bayes_cor (r_cro_mat)[1]

class Bayes:
	def __init__(self, bayesType):
		self.name = "Bayes%i" % bayesType
		self.parallel = True
		self.bayesType = bayesType


	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot)
		global current_bayes_cor, cor1, cor2, cor3
		if self.bayesType == 1:
			current_bayes_cor = cor1
		elif self.bayesType == 2:
			current_bayes_cor = cor2
		elif self.bayesType == 3:
			current_bayes_cor = cor3
		else:
			print "Invalid bayes selection"
			sys.exit()
			#TODO row error instead of print


	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

	# have this call always after __init__ since init initialize the right bayes function
	calculateScore = staticmethod(calculateScoreBayes)

# @ author Florian Goebels
# returns weighted cross correlation similarity score which is based on previously published methods Pierre 2011 et al.
def calculateScore_wcc(a,b):
	global r_wcc
	return r_wcc(robjects.FloatVector(a), robjects.FloatVector(b), 1)[0]

class Wcc:
	calculateScore = staticmethod(calculateScore_wcc)

	def __init__(self):
		self.name="wcc"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True


# @ author Florian Goebels
# returns travor correlation which is pearson correlation plus poisson noise to remove the influence of low counts
def calculateScore_PCCPN(a,b):
	global Poisson_cor_Mat
	index, a = a
	_, b = b
	return Poisson_cor_Mat[index][a][b]

def traver_corr(mat, repeat=1000, norm='columns', verbose=True):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each column to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat, iteration_display, norm):
        if verbose: print iteration_display
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        if norm=='columns':
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
        elif norm=='rows': # seems to make no performance difference 1/25
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
        corr = np.nan_to_num(np.corrcoef(poisson_mat))
        return corr
    avg_result = (reduce(operator.add, (poisson_corr(mat, i, norm=norm) for i in
                                        range(repeat))) / repeat)
    return avg_result

class Poisson:
	calculateScore = staticmethod(calculateScore_PCCPN)

	def __init__(self, repeat=100):
		self.name="poisson-%i" % (repeat)
		self.repeat=repeat
		self.parallel = True

	def init(self, elutionData):
		global Poisson_cor_Mat
		global prot2profile
		fname = "%s.%s" % (elutionData.name, self.name)
		noise_cor_mat = traver_corr(elutionData.elutionMat, self.repeat, 'columns', True)
		Poisson_cor_Mat.append(noise_cor_mat)
		index = len(Poisson_cor_Mat)-1
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = (index, elutionData.getProtIndex(prot))

	def setPoisson_cor_Mat(self, elutionData):
#		Poisson_cor_MatTasks = []
#		global num_cores
#		num_jobs = num_cores
#		for i in range(num_jobs):
#			Poisson_cor_MatTasks.append(tuple([traver_corr, tuple([elutionData.elutionMat, int(self.repeat/num_jobs), 'columns', True])]))
#		p = mp.Pool(num_cores)
#		tmpMats = p.map(no_queue_worker, Poisson_cor_MatTasks)
#		p.close()
#		p.join()
#
		global Poisson_cor_Mat
#		Poisson_cor_Mat = (reduce(operator.add, tmpMats)/self.repeat)
		Poisson_cor_Mat = traver_corr(elutionData.elutionMat, self.repeat, 'columns', True)

	def clear(self):
		global Poisson_cor_Mat
		Poisson_cor_Mat = ""
		global prot2profile
		prot2profile = {}

# @ author Florian Goebels
# returns Mutual Information of two proteins
# mutual information is based on entropy calculation MI(x,y) = H(x) + H(y) - H(x,y) 
# wehre H means entropy
def calculateScore_MI(a, b):
	entropy_a, a_upper, a_lower, _ = a
	entropy_b, b_upper, b_lower, numFracs = b
	joint_probs = np.array(map(lambda x: len(x)/numFracs, [a_upper&b_upper, a_upper&b_lower, a_lower&b_upper, a_lower&b_lower]))
	joint_entropy_a_b = entropy(joint_probs, 2)
	mutual_information =  entropy_a  + entropy_b - joint_entropy_a_b
	return mutual_information

def bin_entropy(p):
	return entropy(np.array([p, 1 - p]))


def entropy(probs, base=0):
	if base == 0: base = len(probs)
	tmp_probs = probs
	tmp_probs[tmp_probs == 0] = 1
	return -sum(probs * map(lambda x: math.log(x, base), tmp_probs))


def getFracs(a, cutoff):
	upper = set([i for i, v in enumerate(a) if v > cutoff])
	lower = set([i for i, v in enumerate(a) if v <= cutoff])
	return (upper, lower)

class MutualInformation():
	calculateScore = staticmethod(calculateScore_MI)

	def __init__(self, minCounts = 2):
		self.name="MI"
		self.minCounts = minCounts
		self.parallel = True
	
	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			profile = elutionData.getElution(prot)
			(profile_upper, profile_lower) = getFracs(profile, self.minCounts)
			numFracs = len(profile)
			prot_entropy = bin_entropy(len(profile_upper) / numFracs)
			prot2profile[fname][prot] = (prot_entropy, profile_upper, profile_lower, numFracs)



	def clear(self):
		global prot2profile
		prot2profile = {}
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
	calculateScore = staticmethod(calculateScore_Jaccard)

        def __init__(self):
                self.name="Jaccard"
		self.non_zero_fracs_for_prot = {}
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = set(np.nonzero(elutionData.getElution(prot))[0])

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True



# @ author Florian Goebels
# return Pearson correlation of two proteins
def calculateScore_Pearson(a,b):
	score = scipy.stats.pearsonr(a, b)[0]
	if math.isnan(score): return 0.0
	return scipy.stats.pearsonr(a, b)[0]

class Pearson:
	calculateScore = staticmethod(calculateScore_Pearson)

	def __init__(self):
		self.name = "Pearson"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot)

	def clear(self):
		global prot2profile
		prot2profile = {}
		return True

# @ author Florian Goebels
# returns Euclidean distance of two proteins
def calculateScore_euclidean(a,b):
	return 1-distance.euclidean(a,b)

class Euclidiean:
	calculateScore = staticmethod(calculateScore_euclidean)

	def __init__(self):
		self.name = "Euclidiean"
		self.parallel = True

	def init(self, elutionData):
		fname = "%s.%s" % (elutionData.name, self.name)
		global prot2profile
		prot2profile[fname] = {}
		for prot in elutionData.prot2Index:
			prot2profile[fname][prot] = elutionData.getElution(prot, normed=True)

	def clear(self):
		global prot2profile
		prot2profile = {}
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
		self.scores = ""
		self.header = ["ProtA","ProtB"]
		self.ppiToIndex = {}

	def addLabels(self, positive, negative):
		self.positive = set([])
		self.negative = set([])

		for ppi in positive:
			if ppi in self.ppiToIndex: self.positive.add(ppi)

		for ppi in negative:
			if ppi in self.ppiToIndex: self.negative.add(ppi)

	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects unto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined tiwth self object
	def merge_singe_ScoreCalc(self, toMerge):

		allPPIs = set([])
		for i in range(self.scores.shape[0]):
			allPPIs.add(self.ppiToIndex[i])

		for j in range(toMerge.scores.shape[0]):
			allPPIs.add(toMerge.ppiToIndex[i])

		numFeature_in_merge = len(toMerge.header)-2
		numFeature_in_self = len(self.header)-2
		new_scores = np.zeros((len(allPPIs), numFeature_in_merge+numFeature_in_self))
		new_ppiToIndex = {}

		k = 0
		for ppi in allPPIs:
			scoresA = [0]*numFeature_in_self
			scoresB = [0]*numFeature_in_merge
			if ppi in self.ppiToIndex:
				scoresA = self.scores[self.ppiToIndex[ppi],:]
			if ppi in toMerge.ppiToIndex:
				scoresB = toMerge.scores[toMerge.ppiToIndex[ppi],:]
			new_score = np.array(list(scoresA) + list(scoresB))
			new_scores[k,:] = new_score
			new_ppiToIndex[k] = ppi
			new_ppiToIndex[ppi] = k
			k+=1
		self.scores = new_scores
		self.ppiToIndex = new_ppiToIndex
		self.header.extend(toMerge.header[2:])


	def getPPIsCoelutionCutoff(self, cutoff = 0.50000, cols = ""):
		if cols == "": cols = range(len(self.header) - 2)
		todel = []
		for i in range(self.scores.shape[0]):
			row = self.scores[i, :]
			hasScore = False
			for score in row:
				if float("{0:.5f}".format(score)) >= cutoff:
					break
			if not hasScore: todel.append(i)
		return todel

	def filter_predictions_coelutionscore(self, cutoff = 0.5, cols = ""):
		todel = self.getPPIsCoelutionCutoff(cutoff, cols)
		self.scores = np.delete(self.scores, todel, 0)
		ppis = np.array(self.ppiToIndex.keys())
		ppis = np.delete(ppis, todel)
		self.ppiToIndex = {}
		for i, ppi in enumerate(ppis):
			self.ppiToIndex[ppi] = i
			self.ppiToIndex[i] = ppi


	def retrieve_scores(self, scores):
		scoreNames = set([])
		score_names_in_header = self.header[2:]
		indicesTokeep = []
		for st in scores: scoreNames.add(st.name)

		for i in range(len(score_names_in_header)):
			head = score_names_in_header[i]
			scorename = head.split(".")[-1]
			if scorename in scoreNames: indicesTokeep.append(i)
		return self.remove_scores(indicesTokeep)

	def remove_scores(self, indices):
		out = CalculateCoElutionScores()
		out.header = ["protA", "protB"]
		for i in indices:
			out.header.append(self.header[i+2])
  		for ppi in self.scores:
			if sum(list(np.array(map(float, self.scores[ppi]))[indices])) !=0:
				out.scores[ppi] = list(np.array(self.scores[ppi])[indices])
		return out

	def rebalance(self, ratio = 5):
		if len(self.positive) * ratio > len(self.negative):
			print "Warning: not enough negative data points in reference to create desired ratio"
		else:
			self.negative = set(random.sample(self.negative, len(self.positive)*ratio))

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
	def filter_interactions_Jaccard(self, eData, ppis):
		print "filtering Jaccard"
		out = set([])
		del_Jaccard_scores = False
		global prot2profile
		if eData.name + ".Jaccard" not in prot2profile:
			this_Jaccard = Jaccard()
			this_Jaccard.init(eData)
			del_Jaccard_scores = True

		for ppi in ppis:
			protA, protB = ppi.split("\t")
			profileA = prot2profile[eData.name + ".Jaccard"][protA]
			profileB = prot2profile[eData.name + ".Jaccard"][protB]
			if not profileA.isdisjoint(profileB): out.add(ppi)

		if del_Jaccard_scores:
			del prot2profile[eData.name + ".Jaccard"]
		return out


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
				protA, protB = sorted([protA, protB])
				allPPIs.add("%s\t%s" % (protA, protB))
		return allPPIs

	# @author: Florian Goebels
	# create co elution table for a given list of scores and ppis
	# @Param:	
	#		scoreTypes	a list of co elutoin score objects
	#		PPIs		a list of ppis for which all scores in scoreTypes schoukd be calculated
	def calculateScores(self, toPred, elutionDatas, scores, outFile):
		task_queue = mp.JoinableQueue()
		out_queue = mp.Queue()
		global num_cores
		num_features = len(scores)*len(elutionDatas)
		for i in range(num_cores):  # remove one core since listener is a sapareted process
			mp.Process(target=getScore, args=(task_queue, out_queue)).start()

		outFH = open(outFile, "w")
		print >> outFH, "\t".join(self.header)
		for ppi in toPred:
			i = 0
			for _ in elutionDatas:
				for score in scores:
					fname = self.header[i+2]
					task = (ppi, i, fname, score.calculateScore)
					task_queue.put(task)
					i += 1
			task_queue.join()
			ppi_scores = [0]*num_features #np.frombuffer(out_array.get_obj())
			for _ in range(num_features):
				score_index, score = out_queue.get()
				ppi_scores[score_index] = score

			print >> outFH, "%s\t%s" % (ppi, "\t".join(map(str, ppi_scores)))

		print "done calcualting co-elution scores"
		for i in range(num_cores):
			task_queue.put('STOP')
		outFH.close()


		"""
		if scoreType.parallel == True:

			task_queue = mp.JoinableQueue()
			out_array = mp.Array(c_float, len(ppiList))
			# starts process which will be feeded by task_queue
			global num_cores
			for i in range(num_cores):  # remove one core since listener is a sapareted process
				mp.Process(target=getScore, args=(task_queue, out_array, scoreType.calculateScore)).start()

			print len(ppiList)
			for i in range(len(ppiList)):
				task_queue.put(i)

			#Wait for all PPIs to be calcualted
			task_queue.join()
			#Wait until all results habe been written to output
			print "done calcualting: " + name
			for i in range(num_cores):
				task_queue.put('STOP')

			for i in range(len(ppiList)):
				self.scores[self.ppiIndex[ppiList[i]]][col_index] = out_array[i]

		else:
			print "Currently only supprot co-elution features than can be calcualted in parallel"
			sys.exit()

		scoreType.clear()
		"""

	def getAllPairs_coelutionDatas(self, elutionDatas):
		allfilteredPPIs = set([])
		for elutionData in elutionDatas:
			print "Filtering: %s" % (elutionData.name)
			scoreCalc = CalculateCoElutionScores(elutionData)
			candidatePPis = scoreCalc.getAllPairs()
			print "Before filtering %i PPIs" % (len(candidatePPis))
			filteredPPIs = scoreCalc.filter_interactions_Jaccard(elutionData, candidatePPis)
			del candidatePPis
			print "After filtering %i PPIs" % (len(filteredPPIs))
			allfilteredPPIs |= filteredPPIs
			del scoreCalc
		print "Num of PPIs across all data stes after filtering %i" % (len(allfilteredPPIs))
		return allfilteredPPIs

	def calculate_coelutionDatas(self, elutionDatas, scores, outDir):
		toPred = self.getAllPairs_coelutionDatas(elutionDatas)
		for eD in elutionDatas:
			for score in scores:
				self.header.append("%s.%s" % (eD.name, score.name))
		self.calculateScores(toPred, elutionDatas, scores, outDir + ".scores.txt")


	# @author: Florian Goebels
	# prints table
	# @Param:
	#		labels print class lable or not
	def toTable(self):
		out = "\t".join(self.header)
		for i in range(self.scores.shape[0]):
			ppi = self.ppiIndex[i]
			scores = "\t".join(map(str, self.scores[i, :]))
			out += "\n%s\t%s" % (ppi, scores)
		return out

	def countLabels(self):
		counts = {}
		for (_, _, label) in self.scores:
			if label not in counts: counts[label] = 0
			counts[label] += 1
		return counts

	def readTable(self, scoreF, cutoff=0.5):
		scoreFH = open(scoreF)
		self.header = scoreFH.readline().rstrip().split("\t")
		i = 0
		self.ppiToIndex = {}
		self.scores = []
		for line in scoreFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			edge = "\t".join(sorted(linesplit[0:2]))
			scores = np.nan_to_num(np.array(map(float, linesplit[2:])))
			if len(np.where(scores > cutoff)[0])<1: continue
			self.ppiToIndex[edge] = i
			self.ppiToIndex[i] = edge
			self.scores.extend(scores)
			i += 1
		scoreFH.close()
		self.scores = np.reshape(self.scores, (i, len(self.header)-2))

	# @author: Florian Goebels
	# prints arff table for weka
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
	def toSklearnData(self, get_preds = True):
		data = []
		targets = []
		ids = []
		for i in range(self.scores.shape[0]):
			ppi = self.ppiToIndex[i]
			label = "?"
			if ppi in self.positive: label = 1
			if ppi in self.negative: label = 0
			scores = self.scores[i,:]
			if label != "?" or get_preds:
				targets.append(label)
				ids.append(ppi)
				data.append(scores)
		return ids, np.array(data), np.array(targets)

# @ author Florian Goebels
# wrapper for machine learning
class CLF_Wrapper:
	# @author: Florian Goebels
	# class initializer, supports both random forest and svm
	# @Param:
	#		data matrix with features where each row is a data point
	#		targets list with class lables
	# 		forest if true ml is random forst, if false ml is svm
	def __init__(self, data, targets, forest=False, folds =10, useFeatureSelection= False):
		self.data = data
		self.folds = folds
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		self.eval_preds = ""
		self.eval_targets = ""
		thisCLF = ""
		if forest:
			print ("using Random forest")
			global num_cores
			thisCLF = RandomForestClassifier(n_estimators=400, n_jobs=num_cores)
		else:	
			print ("Using SVM")
			thisCLF =  svm.SVC(kernel="linear", probability=True)

		if useFeatureSelection:
			self.clf = Pipeline([
				('feature_selection', RFECV(estimator = thisCLF, step = 1, scoring="accuracy")),
				('classification', thisCLF)
			])
		else:
			self.clf = thisCLF

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
	def predict_proba(self, toPred):
		probas = self.clf.predict_proba(toPred)
		return probas[:,1]

	def predict(self, toPred):
		preds = self.clf.predict(toPred)
		return preds

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
	for (name, curve) in curves:
		x, y = curve[0:2]
		if name != "":
			plt.plot(x, y, label=name, color = cols.pop())
		else:
			plt.plot(x, y, color=cols.pop())
	art = []
	lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1),  ncol = 5, fontsize=8)
	art.append(lgd)
	plt.savefig(outF, additional_artists=art, bbox_inches="tight")

# @author Florian Goebels
def predictInteractions(scoreCalc, outDir, useForest):
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(get_preds=False)
	print data_train.shape
	ids_pred, data_pred, targets_pred = scoreCalc.toSklearnData()

	print data_pred.shape
	clf = CLF_Wrapper(data_train, targets_train, forest=useForest, useFeatureSelection=False)
	pred_prob = clf.predict_proba(data_pred)
	pred_class = clf.predict(data_pred)

	outFH = open(outDir + ".pred.txt", "w")
	for i in range(len(ids_pred)):
		if pred_class [i] == 1:
			outFH.write("%s\t%f\n" % ("\t".join(ids_pred[i]), pred_prob[i]))

	pred_prob = clf.predict_proba(data_train)
	pred_class = clf.predict(data_train)
	for i in range(len(ids_train)):
		if pred_class [i] == 1:
			outFH.write("%s\t%f\n" % ("\t".join(ids_pred[i]), pred_prob[i]))

	outFH.close()




def get_eval(scoreCalc, useForest=False, folds=10):
	_, data, targets = scoreCalc.toSklearnData(get_preds=False)
	data = np.array(data)
	print data.shape
	print targets
	clf = CLF_Wrapper(data, targets, forest=useForest, folds=folds)
	eval_scores = clf.getValScores()
	p, r, tr = clf.getPRcurve()
	pr_curve = tuple([r,p, tr])
	fpr, tpr, tr = clf.getROCcurve()
	roc_curve = tuple([fpr, tpr, tr])
	return tuple([eval_scores, pr_curve, roc_curve])

def bench_scores(scoreCalc, outDir, useForest=False, folds=10):
	pr_curves = []
	roc_curves = []
	eval_scores = []
	cutoff_curves = []
	es, pr, roc = get_eval(scoreCalc, useForest, folds)
	pr_curves.append(("Combined", pr))
	roc_curves.append(("Coxmbined", roc))
	eval_scores.append(("Combined", es))
	plotCurves(pr_curves, outDir + ".pr.pdf", "Recall", "Precision")
	plotCurves(roc_curves, outDir + ".roc.pdf", "False Positive rate", "True Positive Rate")
	recall, precision, threshold = pr
	threshold = np.append(threshold, 1)
	cutoff_curves.append(("Precision", (precision, threshold)))
	cutoff_curves.append(("Recall", (recall, threshold)))
	plotCurves(cutoff_curves, outDir + ".cutoff.pdf", "Cutoff", "Evaluation metric score")

	tableFH = open(outDir + ".eval.txt", "w")
	tableFH.write("Name\tPrecision\tRecall\tF-measure\tau_PR\tau_ROC\n")
	for (score_name, evals) in eval_scores:
		tableFH.write("%s\t%s\n" % (score_name, "\t".join(map(str, evals))))
	tableFH.close()

# @author Florian Goebels
def main():
	refF, useForest, elutionFiles, outDir = sys.argv[1:]
	useForest = useForest == "True"

#	scores = [MutualInformation(2), Bayes3(), Euclidiean(), Wcc(), Jaccard(), Poisson(100)]
	scores = [Bayes(3), Wcc(), Poisson(50) , MutualInformation(2), Euclidiean(), Jaccard(), Apex(), Pearson()]
#	scores = [ Wcc(), Jaccard(), Euclidiean(), Bayes(3), MutualInformation(2)]
#	scores = [Poisson(50), Wcc()]

	elutionFH = open(elutionFiles)
	elutionDatas = []
	elutionProts = set([])
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())
		for score in scores:
			score.init(elutionData)

	numFeatures = len(elutionDatas) * len(scores)
#	Calcualte reference scores
	reference = GS.Goldstandard_from_reference_File(refF, found_prots=elutionProts)
	positive = reference.goldstandard_positive
	negative = reference.goldstandard_negative

#	reference_from_human = GS.Goldstandard_from_CORUM("10090", found_prots=elutionProts) # Human: 9606, Worm: 6239,
#	reference_from_mouse = GS.Goldstandard_from_CORUM("10090", found_prots=elutionProts, source_species_regex = "Mouse") # Human: 9606, Worm: 6239,
#	positive =  reference_from_human.goldstandard_positive | reference_from_mouse.goldstandard_positive
#	negative = reference_from_human.goldstandard_negative | reference_from_mouse.goldstandard_negative

	scoreCalc = CalculateCoElutionScores()
	scoreCalc.calculate_coelutionDatas(elutionDatas, scores, outDir)

	global prot2profile
	del prot2profile
	scoreCalc = CalculateCoElutionScores()
	scoreCalc.readTable(outDir + ".scores.txt", 0.5)

#	reference = GS.GS_from_PPIs(refF, found_ppis=set(scoreCalc.ppiToIndex.keys()))
#	positive = reference.goldstandard_positive
#	negative = reference.goldstandard_negative


	scoreCalc.addLabels(positive, negative)


	scoreCalc.rebalance()
#	print "Done"
	print len(scoreCalc.positive)
	print len(scoreCalc.negative)

	bench_scores(scoreCalc, outDir, useForest)
	sys.exit()
	predictInteractions(scoreCalc, outDir, useForest)


if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
