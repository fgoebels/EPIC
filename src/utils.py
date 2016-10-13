#!/usr/bin/env python
from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import copy
import numpy as np
import sys
import os


def getOverlapp(complexesA, complexesB):
	out = 0
	for comp_ID_A in complexesA:
		protsA = complexesA[comp_ID_A]
		matched = False
		for comp_ID_B in complexesB:
			protsB = complexesB[comp_ID_B]
			if (len(protsA & protsB) / min(len(protsA),len(protsB))) > 0.5:
				matched = True
				break

		if matched: out += 1
	return out

def readclusters(cluster_file):
	out = {}
	cluster_FH = open(cluster_file)
	i = 0
	for line in cluster_FH:
		line = line.rstrip()
		prots = set(line.split("\t"))
		if len(prots)<3:continue
		out[i] = prots
		i +=1
	cluster_FH.close()
	return out


def get_cluster_overlapp():

	cluster_dir, outF = sys.argv[1:]

	all_clusters = []

	for file in os.listdir(cluster_dir):
		if file.startswith(".") : continue
		clusters = readclusters(cluster_dir + os.sep + file)
		clustername = file
		all_clusters.append((clustername, clusters))

	out = ""
	for nameA, clustersA, in all_clusters:
		out += nameA
		for nameB, clustersB in all_clusters:
			if nameA == nameB:
				out += "\t%i" % (len(clustersA))
				continue
			out += "\t%i" % (getOverlapp(clustersA, clustersB))
		out += "\n"
	outFH = open(outF, "w")
	print >> outFH, out
	outFH.close()

def cluster_overlapp():

	reference_corum = GS.Goldstandard_from_reference_File("/Users/florian/Desktop/Ce_gs.txt", found_prots="")
	positive_corum = reference_corum.goldstandard_positive
	negative_corum = reference_corum.goldstandard_negative

	reference_go = GS.Goldstandard_from_cluster_File("/Users/florian/workspace/scratch/EPIC_out/clusters/Go_elegans_complex_experiments_mapped.txt", found_prots="")
	positive_go = reference_go.goldstandard_positive
	negative_go = reference_go.goldstandard_negative

	combined = positive_corum&positive_go

	print "\tCORUM\tGO\tCORUM+GO"
	print "positive\t%i\t%i\t%i" % (len(positive_corum - combined), len(positive_go - combined), len(positive_corum&positive_go))
	print "negative\t%i\t%i\t%i" % (len(negative_corum - combined), len(negative_go - combined), len(negative_corum&negative_go))

	scorecalc = CS.CalculateCoElutionScores()
	scorecalc.readTable("/Users/florian/workspace/scratch/EPIC_out/All.scores.sel.txt")
	found_ppis = set(scorecalc.ppiToIndex.keys())
	positive_go     = found_ppis & positive_go
	negative_go     = found_ppis & negative_go
	negative_corum  = found_ppis & negative_corum
	positive_corum  = found_ppis & positive_corum
	combined = positive_corum & positive_go

	print "\tCORUM\tGO\tCORUM+GO"
	print "positive\t%i\t%i\t%i" % (len(positive_corum - combined), len(positive_go - combined), len(positive_corum&positive_go))
	print "negative\t%i\t%i\t%i" % (len(negative_corum - combined), len(negative_go - combined), len(negative_corum&negative_go))

def main():
	get_cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
