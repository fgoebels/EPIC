from __future__ import division

import mmap
import numpy as np
import scipy.stats
import sys
import math
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from scipy.spatial import distance
import operator
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score, precision_recall_fscore_support
from sklearn.model_selection import  cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFECV
import matplotlib.pyplot as plt
import inspect
import os
from sklearn import svm
from sklearn import metrics
import random
import multiprocessing as mp
import GoldStandard as GS

# The following imported libraries are for intergrating GeneMANIA data as Functional Evidence
import glob
import sys
from collections import defaultdict
from bs4 import BeautifulSoup
import urllib
import urllib2
import os 


# GLobal static objects used across the files
# Load R packages globally in order to avoid loading each package again in each thread

# Load script for calculating Bayes correlation globally
r=robjects.r
print os.path.realpath(__file__)
r.source(os.path.realpath(__file__).rsplit(os.sep,1)[0] + os.sep + "Bayes_Corr.R")

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

def lineCount(filename):
	f = open(filename, "r+")
	buf = mmap.mmap(f.fileno(), 0)
	lines = 0
	readline = buf.readline
	while readline():
		lines += 1
	return lines


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
			prot2profile[fname][prot] = map( math.log, elutionData.getElution(prot)+1)
		global current_bayes_cor, cor1, cor2, cor3
		if self.bayesType == 1:
			current_bayes_cor = cor1
		elif self.bayesType == 2:
			current_bayes_cor = cor2
		elif self.bayesType == 3:
			current_bayes_cor = cor3
		else:
			print("Invalid bayes selection")
			sys.exit()
			#TODO throw error instead of print


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
		if verbose: print(iteration_display)
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
	
# @ author Lucas Ming Hu
# This is a helper class for getting GeneMANIA functional evidence for a given ElutionData object
class Genemania:
	
	def __init__(self, taxID):
		self.taxoID = taxID
		# Get all Genemania files
		self.catchFile()
		# all functional evidence codes in GeneMANIA, excluding "Physical" and "complexes" and "Predicted" to eliminate circularity
		self.functionalEvidences = ['Co-expression', 'Genetic', 'Other', 'Shared']
		# loads all of Worm Gene
		self.load_genemania()

#	def get_score(self, ppi):

	# @auothor Lucas Ming Hu
	# the catchFile function can help to download files from GeneMANIA website automatically.	
	def catchFile(self): 
		taxoIDspeciesDic = {'3702':'Arabidopsis_thaliana', '6239':'Caenorhabditis_elegans', '7955':'Danio_rerio', 
		                    '7227':'Drosophila_melanogaster','562':'Escherichia_coli','9606':'Homo_sapiens',
		                    '10090':'Mus_musculus','10116':'Rattus_norvegicus','4932':'Saccharomyces_cerevisiae'} 
		
		if self.taxoID not in taxoIDspeciesDic:
			return None #TODO throw illegal argument exception
	
		urlbase = 'http://genemania.org/data/current'
		speciesURL = os.path.join(urlbase, taxoIDspeciesDic[self.taxoID])
		r = urllib.urlopen(speciesURL).read()
		soup = BeautifulSoup(r)
    
		table = soup.find('table')
    
		allcell = []
		for row in table.find_all('tr'):
			for col in row.find_all('td'):
				allcell.append(col.getText())
    
		#filtering 
		self.files = []
		for c in allcell:
			if '.txt' in c:
				self.files.append(os.path.join(speciesURL,c))

	# @author: Lucas Ming Hu        
	# a helper function to get the average of the GeneMANIA scores 
	# for each line of evidence
	def average(self, secondaryEvidenceDic):
		resultDict = defaultdict(float)		
		for key in secondaryEvidenceDic:
			resultDict[key] = sum(secondaryEvidenceDic[key]) * 1.0 / len(secondaryEvidenceDic[key])
		return resultDict
	
	
	# returns Functional anotation scores as a CalculateCoElutionScores Object
	def load_genemania(self):
		self.scores = {}
		self.header = []
		# read online database - by species taxoID

		for i, f_evidence in enumerate(self.functionalEvidences):
			this_evidence_scores = {}
			self.header.append("GeneMania_%s" % f_evidence)
			for fp in self.files:
				filename = str(fp.split('/')[-1])
				if filename.startswith(f_evidence):
					print "Processing: %s" % (filename)
					fh = urllib2.urlopen(fp)
					fh.readline()
					for line in fh:
						proteinA, proteinB, score = line.split()
						edge = "\t".join(sorted([proteinA, proteinB]))
						score = float(score)
						if edge not in this_evidence_scores: this_evidence_scores[edge] = np.array([0, 0])
						this_evidence_scores[edge][0] += score
						this_evidence_scores[edge][1] += 1
					fh.close()

		for edge in this_evidence_scores:
			score, counts = this_evidence_scores[edge]
			avg_score = score/counts
			if edge not in self.scores: self.scores[edge] = [0]*len(self.functionalEvidences)
			self.scores[edge][i] = avg_score


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
		self.scores = np.array([])
		self.header = ["ProtA","ProtB"]
		self.ppiToIndex = {}
		self.IndexToPpi = {}
		self.positive = set([])
		self.negative = set([])

	def addLabels(self, positive, negative):
		self.positive = positive
		self.negative = negative

#		for ppi in positive:
#			if ppi in self.ppiToIndex: self.positive.add(ppi)

#		for ppi in negative:
#			if ppi in self.ppiToIndex: self.negative.add(ppi)


	# @author: Florian Goebels
	# this method combines who CalculateCoElutionScores objects unto one by comping the toMerge object into the self object
	# @Param:
	#		CalculateCoElutionScores toMerge a second CalculateCoElutionScores which should be combined tiwth self object
	def merge_singe_ScoreCalc(self, toMerge):
		allPPIs = set(self.ppiToIndex.keys()) | set(toMerge.ppiToIndex.keys())
		numFeature_in_merge = len(toMerge.header)-2
		numFeature_in_self = len(self.header)-2
		new_scores = np.zeros((len(allPPIs), numFeature_in_merge+numFeature_in_self))
		new_ppiToIndex = {}
		new_IndexToPpi = {}
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
			new_IndexToPpi[k] = ppi
			new_ppiToIndex[ppi] = k
			k+=1
		self.scores = new_scores
		self.ppiToIndex = new_ppiToIndex
		self.IndexToPpi = new_IndexToPpi
		self.header.extend(toMerge.header[2:])

	def filter_coelutionscore(self, score_cutoff = 0.5):
		valid_rows = []
		dim = self.scores.shape
		for i in range(dim[0]):
			for j in range(dim[1]):
				if self.scores[i][j] > 0.5:
					valid_rows.append(i)
					break
#		valid_rows = list(set(np.where(self.scores>score_cutoff)[0]))
		self.scores = self.scores[valid_rows,:]
		new_IndexToPpi = {}
		new_ppiToIndex = {}
		for i in range(len(valid_rows)):
			old_index = valid_rows[i]
			ppi = self.IndexToPpi[old_index]
			new_ppiToIndex[ppi] = i
			new_IndexToPpi[i] = ppi
		self.IndexToPpi = new_IndexToPpi
		self.ppiToIndex = new_ppiToIndex

	def retrieve_scores(self, scores):
		scoreNames = set([])
		for st in scores: scoreNames.add(st.name)
		to_keep_header = [0,1]
		to_keep_score = []
		for i in range(2,len(self.header)):
			colname = self.header[i]
			scorename = colname.split(".")[-1]
			if scorename in scoreNames:
				to_keep_header.append(i)
				to_keep_score.append(i-2)

		self.header = list(np.array(self.header)[to_keep_header])
		self.scores = self.scores[:,to_keep_score]
		self.filter_coelutionscore()

	def rebalance(self, ratio = 5):
		if len(self.positive) * ratio > len(self.negative):
			print("Warning: not enough negative data points in reference to create desired ratio")
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
		print("filtering Jaccard")
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
	def calculateScores(self, toPred, elutionDatas, scores, outFile, num_cores, verbose = True):
		task_queue = mp.JoinableQueue()
		out_queue = mp.Queue()
		self.scores = []
		self.ppiToIndex = {}
		gs = self.positive | self.negative
		print len(gs)
		num_features = len(scores)*len(elutionDatas)
		self.scores = np.zeros((len(gs), num_features))
		for i in range(num_cores):  # remove one core since listener is a sapareted process
			mp.Process(target=getScore, args=(task_queue, out_queue)).start()
		outFH = open(outFile, "w")
		print >> outFH, "\t".join(self.header)
		k = 0
		ppi_index = 0
		write_buffer = ""
		for ppi in toPred:
			k += 1
			if k % 100000 == 0:
				if verbose: print(k)
				outFH.write(write_buffer)
				outFH.flush()
				write_buffer = ""
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
			ppi_scores = np.nan_to_num(np.array(ppi_scores))
			if len(list(set(np.where(ppi_scores > 0.5)[0]))) > 0:
				write_buffer +=  "%s\t%s\n" % (ppi, "\t".join(map(str, ppi_scores)))
				if ppi in gs:
					self.ppiToIndex[ppi] = ppi_index
					self.IndexToPpi[ppi_index] = ppi
					self.scores[ppi_index,:] = ppi_scores
					ppi_index += 1
		print >> outFH, write_buffer
		print("done calcualting co-elution scores")
		print ppi_index
		self.scores = self.scores[0:ppi_index,:]
		for i in range(num_cores):
			task_queue.put('STOP')
		outFH.close()

	def getAllPairs_coelutionDatas(self, elutionDatas):
		allfilteredPPIs = set([])
		for elutionData in elutionDatas:
			print("Filtering: %s" % (elutionData.name))
			scoreCalc = CalculateCoElutionScores(elutionData)
			candidatePPis = scoreCalc.getAllPairs()
			print("Before filtering %i PPIs" % (len(candidatePPis)))
			filteredPPIs = scoreCalc.filter_interactions_Jaccard(elutionData, candidatePPis)
			del candidatePPis
			print("After filtering %i PPIs" % (len(filteredPPIs)))
			allfilteredPPIs |= filteredPPIs
			del scoreCalc
		print("Num of PPIs across all data stes after filtering %i" % (len(allfilteredPPIs)))
		return allfilteredPPIs

	def calculate_coelutionDatas(self, elutionDatas, scores, outDir, num_cores, toPred=""):
		if toPred == "":
			toPred = self.getAllPairs_coelutionDatas(elutionDatas)
		for eD in elutionDatas:
			for score in scores:
				self.header.append("%s.%s" % (eD.name, score.name))
		self.calculateScores(toPred, elutionDatas, scores, outDir + ".scores.txt", num_cores)

	# @author: Florian Goebels
	# prints table
	# @Param:
	#		labels print class lable or not
	def toTable(self, fh="", labels= True):
		valid_ppis = set(self.ppiToIndex.keys())
		if labels: valid_ppis = valid_ppis & (self.positive | self.negative)
		out = ""
		if fh!="":
			print >> fh, "\t".join(self.header)
		else:
			out = "\t".join(self.header)
		for i in range(self.scores.shape[0]):
			ppi = self.IndexToPpi[i]
			if ppi not in valid_ppis: continue
			scores = "\t".join(map(str, self.scores[i, :]))
			line = "%s\t%s" % (ppi, scores)
			if fh != "":
				print >> fh, line
			else:
				out += "\n" + line
		if fh !="":fh.flush()
		return out

	def countLabels(self):
		counts = {}
		for (_, _, label) in self.scores:
			if label not in counts: counts[label] = 0
			counts[label] += 1
		return counts

	def readTable(self, scoreF):
		scoreFH = open(scoreF)
		self.header = scoreFH.readline().rstrip().split("\t")
		gs = self.positive | self.negative
		self.scores = np.zeros((len(self.positive)+len(self.negative), len(self.header)-2))
		i = 0
		self.ppiToIndex = {}
		for line in scoreFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			edge = "\t".join(sorted(linesplit[0:2]))
			if edge not in gs: continue
			edge_scores = np.nan_to_num(np.array(map(float, linesplit[2:])))
			self.scores[i,:] = edge_scores
			self.IndexToPpi[i] = edge
			self.ppiToIndex[edge] = i
			i += 1
		scoreFH.close()
		self.scores = self.scores[0:i, :]

	# @author: Florian Goebels
	# return stored co elution scores for a given data set in form that is required for sklearn to learn and predict interaction
	def toSklearnData(self, get_preds = True):
		ids = []
		targets = []
		used_indeces = []
		for i in range(self.scores.shape[0]):
			ppi = self.IndexToPpi[i]
			label = "?"
			if ppi in self.positive: label = 1
			if ppi in self.negative: label = 0
			if label != "?" or get_preds:
				targets.append(label)
				ids.append(ppi)
				used_indeces.append(i)
		if get_preds == True:
			return ids, self.scores, np.array(targets)
		else:
			return ids, self.scores[used_indeces, :], np.array(targets)

# @ author Florian Goebels
# wrapper for machine learning
class CLF_Wrapper:
	# @author: Florian Goebels
	# class initializer, supports both random forest and svm
	# @Param:
	#		data matrix with features where each row is a data point
	#		targets list with class lables
	# 		forest if true ml is random forst, if false ml is svm
	def __init__(self, data, targets, num_cores, forest=False, folds =10, useFeatureSelection= False):
		self.data = data
		self.folds = folds
		self.targets = targets #label_binarize(targets, classes=[0, 1])
		self.eval_preds = ""
		self.eval_targets = ""
		self.num_cores = num_cores
		thisCLF = ""
		if forest:
#			print("using Random forest")
			thisCLF = RandomForestClassifier(n_estimators=400, n_jobs=self.num_cores)
		else:	
#			print("Using SVM")
			thisCLF =  svm.SVC(kernel="linear", probability=True)
			if useFeatureSelection:
				self.clf = Pipeline([
					('feature_selection', RFECV(estimator=thisCLF, step=1, scoring="accuracy")),
					('classification', thisCLF)
				])
		self.clf = thisCLF
		self.clf.fit(self.data, self.targets)

	# @author Florian Goebels
	def geteval(self):
		all_probas = []
		all_targets = []
		skf = StratifiedKFold(self.folds)
		for train, test in skf.split(self.data, self.targets): #Depricated code StratifiedKFold(self.targets, self.folds):
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
		skf = StratifiedKFold(self.folds)
		folds = skf.split(self.data, self.targets) # Depricated code StratifiedKFold(self.targets, self.folds)
		preds = cross_val_predict(self.clf, self.data, self.targets, cv=folds, n_jobs=self.num_cores)
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
def predictInteractions(scoreCalc, outDir, useForest, num_cores, verbose= False):
	All_score_FH = open(outDir + ".scores.txt")


	ids_train, data_train, targets_train = scoreCalc.toSklearnData(get_preds=False)
	clf = CLF_Wrapper(data_train, targets_train, num_cores=num_cores, forest=useForest, useFeatureSelection=False)
	print len(ids_train)
	print data_train.shape
	def getPredictions(scores, edges, clf):
		out = []
		pred_prob = clf.predict_proba(scores)
		pred_class = clf.predict(scores)
		for i, prediction in enumerate(pred_class):
			if prediction == 1:
				out.append("%s\t%f" % (edges[i], pred_prob[i]))
		return out
	out = []
	tmpscores = np.zeros((100000, data_train.shape[1]))
	edges = [""]*100000
	All_score_FH.readline()
	k = 0
	for line in All_score_FH:
		if k % 100000==0:
			out.extend(getPredictions(tmpscores, edges, clf))
			tmpscores = np.zeros((100000, data_train.shape[1]))
			k = 0
		line = line.rstrip()
		if line =="":continue
		linesplit = line.split("\t")
		edge = "\t".join(sorted(linesplit[0:2]))
		edge_scores = np.nan_to_num(np.array(map(float, np.array(linesplit[2:]), ))).reshape(1, -1)
		edges[k] = edge
		tmpscores[k,:] = edge_scores
		k += 1
	out.extend(getPredictions(tmpscores[0:k,:], edges[0:k], clf))
	All_score_FH.close()

	outFH = open(outDir + ".pred.txt", "w")
	outFH.write("\n".join(out))
	outFH.close()
	return outDir + ".pred.txt", len(out)


def load_data(data_dir, scores):
	paths = [os.path.join(data_dir,fn) for fn in next(os.walk(data_dir))[2]]
	elutionDatas = []
	elutionProts = set([])
	for elutionFile in paths:
		if elutionFile.rsplit(os.sep, 1)[-1].startswith("."): continue
		elutionFile = elutionFile.rstrip()
		elutionData = ElutionData(elutionFile)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())
		for score in scores:
			score.init(elutionData)
	return elutionProts, elutionDatas

def create_goldstandard(target_taxid, valprots):
	def create_gs_set(cluster_obj, target_taxid, name, valprots):
		gs =  GS.Goldstandard_from_Complexes(name)
		gs.make_reference_data(cluster_obj, target_taxid, found_prots=valprots)
		return gs

	# Create gold standard from CORUM
	corum_gs = create_gs_set(GS.CORUM(), target_taxid, "CORUM", valprots)

	# Create gold standard from GO
	go_gs = create_gs_set(GS.QuickGO(target_taxid), target_taxid, "GO", valprots)

	#Create gold stnadard from Intact
	intact_gs = create_gs_set(GS.Intact_clusters(), target_taxid, "Intact", valprots)

	corum_p, corum_n = corum_gs.get_goldstandard()
	intact_p, intact_n = intact_gs.get_goldstandard()
	go_p, go_n = go_gs.get_goldstandard()

	all_p = corum_p | intact_p | go_p
	all_n = corum_n | intact_n | go_n

	training_p = corum_p | intact_p - go_p
	training_n = corum_n | intact_n - go_n

	go_complexes = go_gs.complexes
	corum_complexes = corum_gs.complexes

	return training_p, training_n, all_p, all_n, go_complexes, corum_complexes



def calculate_features():
	return None

def make_classification():
	return None

def main():
	out_dir = "/Users/florian/workspace/scratch/EPIC_out/test_output_dir/test"
	selection = [True, False, True, False, False, False, False, False ]
	use_rf = True
	input_dir = "/Users/florian/workspace/scratch/EPIC_out/test_input_dir"
	num_cores = 4
	target_taxid = "6239"
	all_scores = [Pearson(), Jaccard(), Apex(), MutualInformation(2), Euclidiean(), Wcc(), Bayes(3), Poisson(5)]
	this_scores = []
	for i, selection in enumerate(selection):
		if selection: this_scores.append(all_scores[i])
	foundprots, elution_datas = load_data(input_dir, this_scores)
	evals = create_goldstandard(target_taxid, foundprots)
	training_p, training_n, all_p, all_n, go_complexes, corum_complexes = evals
	scoreCalc = CalculateCoElutionScores()
	scoreCalc.addLabels(all_p, all_n)
	scoreCalc.readTable("/Users/florian/workspace/scratch/EPIC_out/test_output_dir/test.scores.txt")
	#scoreCalc.calculate_coelutionDatas(elution_datas, this_scores, out_dir, num_cores)
	scoreCalc.addLabels(training_p, training_n)
	print "doing benchmark"
	bench_scores(scoreCalc, out_dir, num_cores, useForest=use_rf)
	scoreCalc.addLabels(all_p, all_n)
	print "doing prediction"
	predictInteractions(scoreCalc, out_dir, use_rf, num_cores)
	clustering_evaluation(evals, scoreCalc, out_dir, ",".join([score.name for score in this_scores]), num_cores, use_rf)

def clustering_evaluation(evals, scoreCalc, outDir, feature_combination, number_of_cores, use_random_forest):
	training_p, training_n, all_p, all_n, go_complexes, corum_complexes = evals
	scoreCalc.addLabels(training_p, training_n)
	_, data, targets = scoreCalc.toSklearnData(get_preds=False)
	print data.shape
	num_training_ppi = data.shape[0]
	data = np.array(data)
	clf = CLF_Wrapper(data, targets, num_cores=number_of_cores, forest=use_random_forest, folds=3,
						 useFeatureSelection=False)
	eval_scores = clf.getValScores()
	predF = "%s.pred.txt" % (outDir)
	predicted_ppis = lineCount(predF)
	pred_clusters = GS.Clusters(need_to_be_mapped=False)
	pred_clusters.read_file("%s.clust.txt" % (outDir))
	pred_clusters.filter_complexes()
	pred_clusters = pred_clusters
	corum_scores = "\t".join(map(str, pred_clusters.clus_eval(corum_complexes)))
	go_scores = "\t".join(map(str, pred_clusters.clus_eval(go_complexes)))
	line = "%s\t%i\t%s\t%i\t%i\t%s\t%s" % (
		feature_combination, num_training_ppi, "\t".join(map(str, eval_scores)), predicted_ppis,
		len(pred_clusters.complexes), corum_scores, go_scores)
	linesplit = line.split("\t")
	for i, cat in enumerate(["Features", "PPi in training set", "Precision", "Recall", "F-measure", "auPR", "auROC", "Num predicted PPIs", "Num predicted clusters", "CORUM mmr", "CORUM overlapp", "CORUM simcoe", "CORUM mean_simcoe_overlap", "CORUM sensetivity", "CORUM accuracy", "CORUM sep", "GO mmr", "GO overlapp", "GO simcoe", "GO mean_simcoe_overlap", "GO sensetivity", "GO accuracy", "GO sep"]):
		print "%s\t\t%s" % (cat, linesplit[i])

	outFH = open("%s.eval.txt" % (outDir), "w")
	print >> outFH, line
	outFH.close()


def get_eval(scoreCalc, num_cores, useForest=False, folds=3):
	_, data, targets = scoreCalc.toSklearnData(get_preds=False)
	print data.shape
	data = np.array(data)
	clf = CLF_Wrapper(data, targets, num_cores=num_cores, forest=useForest, folds=folds,
					  useFeatureSelection=False)
	p, r, tr = clf.getPRcurve()
	pr_curve = tuple([r,p, tr])
	fpr, tpr, tr = clf.getROCcurve()
	roc_curve = tuple([fpr, tpr, tr])
	return tuple([pr_curve, roc_curve])

def bench_scores(scoreCalc, outDir, num_cores, useForest=False, folds=3):
	pr_curves = []
	roc_curves = []
	cutoff_curves = []
	pr, roc = get_eval(scoreCalc, num_cores, useForest, folds)
	pr_curves.append(("Combined", pr))
	roc_curves.append(("Coxmbined", roc))
	plotCurves(pr_curves, outDir + ".pr.png", "Recall", "Precision")
	plotCurves(roc_curves, outDir + ".roc.png", "False Positive rate", "True Positive Rate")
	recall, precision, threshold = pr
	threshold = np.append(threshold, 1)
	cutoff_curves.append(("Precision", (precision, threshold)))
	cutoff_curves.append(("Recall", (recall, threshold)))
	plotCurves(cutoff_curves, outDir + ".cutoff.png", "Cutoff", "Evaluation metric score")


if __name__ == "__main__":
		try:
				main()
		except KeyboardInterrupt:
				pass
