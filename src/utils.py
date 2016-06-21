#!/usr/bin/env python
from __future__ import division
import sys
import os
import re
import subprocess
import gzip
import math
import numpy as np


def getSubset(tmpData, keys):
        (a,b) = (set(tmpData[asKey(keys[0])]), set(tmpData[asKey(keys[1])]))
        return mkGroups(a, b)

def mkGroups(a, b):
        both = a & b
        total = a | b
        return (a - both, b - both, both, total)

def makeRtable(phosphoDataF, useNA = True, fastaF = ""):
	out = "PoshpoSiteID"

	if fastaF != "":
		fastaData = readFastaFile(fastaF)
		out += "\tProtSize"
	
	phosData = readData(phosphoDataF, np.array([4,5,3,2]), np.array([0,6,8,9]))


        counts = {}
        for experiment in phosData:
                for (protID, start, ascore, copyNmbr) in phosData[experiment]:
                        prot = "_".join([protID, start])
			if fastaF != "": prot += "\t%i" % (len(fastaData[protID]))
                        if prot not in counts: counts[prot] = {}
                        if experiment not in counts[prot]: counts[prot][experiment] = 0
                        counts[prot][experiment] += int(copyNmbr)


        allExperiments = sorted(phosData.keys())

	

        for experiment in allExperiments:
                expname = "_".join(experiment)
                out += "\t%s" % (expname)
        out += "\n"

        for prot in counts:
                prot
                thisCounts = ""
                for experiment in allExperiments:
                        if experiment not in counts[prot]:
                                if useNA: 
					thisCounts += "\tNA"
				else:
					 thisCounts += "\t0"
                        else:
                                thisCounts += "\t%i" % (counts[prot][experiment])
                out += prot + thisCounts + "\n"

	return out

def getCollumnAsArray(index, matrix):
	out = set([])
	for row in matrix:
		out.add(row[index])
	return list(out)

def getOverlapps(thisSet, dictWithSets, keys = ""):
	out = []
	if keys == "": keys = dictWithSets.keys()
	for key in keys:
		out.append(tuple([key, thisSet & set(dictWithSets[key])]))
	return out

def asKey(key):
	return tuple([key])

def containsReturnResultAsString(key, dic, trueVal = "Y", falseVal = "N"):
	if key in dic:
		return trueVal
	else:
		return falseVal

def con(array, delim = "\t"):
	return delim.join(map(str, array))

def readData(dataF, keyColNumbers, valueRowNumbers, primKeyMap = "", header = True, sep="\t"):
	out = {}
	dataFH = open(dataF)
	if header: dataFH.readline() 
	for line in dataFH:
		line = line.rstrip()
		lineSplit = np.array(line.split(sep))
#		if len(lineSplit)<2: continue
		key = tuple(lineSplit[keyColNumbers])
		primkey = tuple([key[0]])
		if primKeyMap != "" and primkey in primKeyMap:
			if len(primKeyMap[primkey]) ==1:
				newPrimKey = primKeyMap[primkey][0][0]
				key = list(key[1:])
				key.insert(0, newPrimKey)
				key = tuple(key)
#			else:
#				print "no uniq map for " + primkey[0]
#		if primKeyMap != "" and primkey not in primKeyMap:
#			print "no mapping found for " + primkey[0]
		value = tuple(lineSplit[valueRowNumbers])
		if key not in out: out[key] = []
		if value not in out[key]: out[key].append(value)
	dataFH.close()
	return out

def readFastaFile(fastaF):
	out = {}
	seq = ""
	thisID = ""
	inF = open(fastaF)
	for line in inF:
		line = line.rstrip()
		matchedID = re.match(">\w\w\|(.*?)\|(.*?) ", line)
		if(matchedID):
			if(thisID == ""):
				thisID =  matchedID.group(1)
				continue
			else:
				out[thisID] = seq
				thisID =  matchedID.group(1)
				seq = ""
				continue
		seq += line
	# flush last entry
	out[thisID] = seq
	inF.close()
	return out



if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
