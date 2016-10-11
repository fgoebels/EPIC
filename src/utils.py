#!/usr/bin/env python
from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import copy
import numpy as np
import sys
import os

from itertools import chain, combinations

def powerset(iterable):
  xs = list(iterable)
  # note we return an iterator rather than a list
  return chain.from_iterable( combinations(xs,n) for n in range(len(xs)+1) )


def benchmark():
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(50), CS.Pearson(), CS.Apex()]
	all_score_names = set([st.name for st in scores])
	useForest = True
	reference = GS.Goldstandard_from_reference_File("/Users/florian/Desktop/Ce_gs.txt", found_prots="")
	positive = reference.goldstandard_positive
	negative = reference.goldstandard_negative

	reference_go = GS.Goldstandard_from_cluster_File("/Users/florian/workspace/scratch/EPIC_out/clusters/Go_elegans_complex_experiments_mapped.txt", found_prots="")
	positive_go = reference_go.goldstandard_positive
	negative_go = reference_go.goldstandard_negative

	positive = positive - positive_go
	negative = negative - negative_go

	go_clusters = readclusters("/Users/florian/workspace/scratch/EPIC_out/clusters/Go_elegans_complex_experiments_mapped.txt")
	corum_cluster = readclusters("/Users/florian/workspace/scratch/EPIC_out/clusters/CORUM_elegans_complexes.txt")

	scorecalc = CS.CalculateCoElutionScores()
	scorecalc.readTable("/Users/florian/workspace/scratch/EPIC_out/All.scores.sel.txt")
	scorecalc.addLabels(positive, negative)

	outFH = open("/Users/florian/workspace/scratch/EPIC_out/FS.txt", "w")
	print >> outFH, "%s\tTraining size\tPrecision\tRecall\tF-measure\tauPR\tauROC\tNum pred PPIs\tNum pred Clusters\tNum clusters overlapp with CORUM\tNum clusters overlapp with GO" % (" ".join(all_score_names))
	for score_subset in ((CS.Pearson(),), (CS.Poisson(50), CS.Wcc()), (scores)):#list(powerset(scores)):
		if len(score_subset) ==0: continue
		this_sc = copy.deepcopy(scorecalc)
		this_sc.retrieve_scores(score_subset)
		this_score_names = set([st.name for st in score_subset])
		fs = []
		for score_name in all_score_names:
			if score_name in this_score_names:
				fs.append("1")
			else:
				fs.append("0")
		fs = " ".join(fs)
		this_sc.rebalance()
		_, data, targets = this_sc.toSklearnData(get_preds=False)
		num_training_ppi =  data.shape[0]
		data = np.array(data)
		if score_subset == (CS.Apex(),):
			line = "%s\t%i\t%s" % (fs, num_training_ppi, "\t".join(["0"]*5))
			print line
			print >> outFH, line
		clf = CS.CLF_Wrapper(data, targets, forest=useForest, folds=10)
		eval_scores = clf.getValScores()
		predF, predicted_ppis = CS.predictInteractions(this_sc,"/Users/florian/workspace/scratch/EPIC_out/All.scores.eval" , useForest)
		out_dir = "/Users/florian/workspace/scratch/EPIC_out/All"
		clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, out_dir)
		os.system(clustering_CMD)
		pred_clusters = readclusters("%s.clust.txt" % (out_dir))
		matched_corum_clusters = getOverlapp(pred_clusters, corum_cluster)
		matched_go_clusters = getOverlapp(pred_clusters, go_clusters)
		line = "%s\t%i\t%s\t%i\t%i\t%i" % (fs, num_training_ppi, "\t".join(map(str, eval_scores)), predicted_ppis, len(pred_clusters) , matched_corum_clusters, matched_go_clusters)
		print line
		print >> outFH, line
	outFH.close()


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
	benchmark()
	#cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
