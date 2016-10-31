from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import numpy as np
import sys
import os

from itertools import chain, combinations

def powerset(iterable):
  xs = list(iterable)
  # note we return an iterator rather than a list
  return chain.from_iterable( combinations(xs,n) for n in range(len(xs)+1) )


def get_cols(header, feature_scores):
	scoreNames = set([])
	for st in feature_scores: scoreNames.add(st.name)
	to_keep_header = [0, 1]
	to_keep_score = []
	for i in range(2, len(header)):
		colname = header[i]
		scorename = colname.split(".")[-1]
		if scorename in scoreNames:
			to_keep_header.append(i)
			to_keep_score.append(i - 2)
	new_header = list(np.array(header)[to_keep_header])
	return new_header, to_keep_score

def readTable(feature_scores, scoreF, score_cutoff = 0.5):
	out = CS.CalculateCoElutionScores()
	num_ppis_infile = CS.lineCount(scoreF) - 1
	with open(scoreF) as scoreFH:
		header = scoreFH.readline().rstrip().split("\t")
		new_header, scores_to_keep = get_cols(header, feature_scores)
		out.scores = np.zeros((num_ppis_infile, len(new_header) - 2))
		i = 0
		out.ppiToIndex = {}
		for line in scoreFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			edge = "\t".join(sorted(linesplit[0:2]))
			edge_scores = np.nan_to_num(np.array(map(float, np.array(linesplit[2:])[scores_to_keep],)))
			if len(list(set(np.where(edge_scores > score_cutoff)[0])))>0:
				out.scores[i, :] = edge_scores
				out.IndexToPpi[i] = edge
				out.ppiToIndex[edge] = i
				i += 1
	out.scores = out.scores[:i,:] # remove empty zero rows at the end of the matrix
	scoreFH.close()
	return out


# self.scores = np.reshape(self.scores, (i, len(self.header)-2))

def predictInteractions(All_score_F, train_scoreCalc, scores, useForest = True, num_cores =4, score_cutoff=0.5):
	out = []
	All_score_FH = open(All_score_F)
	header = All_score_FH.readline().rstrip().split("\t")
	new_header, scores_to_keep = get_cols(header, scores)

	ids_train, data_train, targets_train = train_scoreCalc.toSklearnData(get_preds=False)
	clf = CS.CLF_Wrapper(data_train, targets_train, num_cores=num_cores, forest=useForest, useFeatureSelection=False)

	for line in All_score_FH:
		line = line.rstrip()
		linesplit = line.split("\t")
		edge = "\t".join(sorted(linesplit[0:2]))
		edge_scores = np.nan_to_num(np.array(map(float, np.array(linesplit[2:])[scores_to_keep], ))).reshape(1, -1)
		if len(list(set(np.where(edge_scores > score_cutoff)[0]))) > 0:
			pred_prob = clf.predict_proba(edge_scores)
			pred_class = clf.predict(edge_scores)
			if pred_class == 1:
				out.append("%s\t%f\n" % (edge, pred_prob))
	All_score_FH.close()
	return out

def benchmark():
	feature_combination, use_random_forest, number_of_cores, refF, train_scoreF, all_scoreF, go_complexesF, corum_complexesF, outDir = sys.argv[1:]
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(50), CS.Pearson(), CS.Apex()]
	use_random_forest = use_random_forest == True
	number_of_cores = int(number_of_cores)
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])

	print this_scores

	scorecalc_train = readTable(this_scores, train_scoreF)

	print scorecalc_train.scores.shape

	reference = GS.Goldstandard_from_reference_File(refF, found_prots="")
	positive = reference.goldstandard_positive
	negative = reference.goldstandard_negative

	reference_go = GS.Goldstandard_from_cluster_File(go_complexesF, found_prots="")
	positive_go = reference_go.goldstandard_positive
	negative_go = reference_go.goldstandard_negative

	positive = positive - positive_go
	negative = negative - negative_go

	go_clusters = readclusters(go_complexesF)
	corum_cluster = readclusters(corum_complexesF)

	scorecalc_train.addLabels(positive, negative)
	scorecalc_train.rebalance()

	_, data, targets = scorecalc_train.toSklearnData(get_preds=False)
	num_training_ppi =  data.shape[0]
	data = np.array(data)
	clf = CS.CLF_Wrapper(data, targets, num_cores=number_of_cores, forest=use_random_forest, folds=2, useFeatureSelection=False)
	eval_scores = clf.getValScores()
	print "Done doing ML"
	predicted_ppis = predictInteractions(all_scoreF, scorecalc_train, this_scores, useForest = use_random_forest, num_cores = number_of_cores)
	predF = "%s.pred.txt" % (outDir)
	predFH = open(predF, "w")
	predFH.write("\n".join(predicted_ppis))
	predFH.close()
	print "Done with prediction"
	clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s%s.clust.txt" % (predF, outDir, feature_combination)
	os.system(clustering_CMD)
	pred_clusters = GS.Clusters(need_to_be_mapped=False)
	pred_clusters.read_file("%s.clust.txt" % (outDir))
	pred_clusters.filter_complexes()
	pred_clusters.merge_complexes()
	pred_clusters = pred_clusters.complexes
	matched_corum_clusters = getOverlapp(pred_clusters, corum_cluster)
	matched_go_clusters = getOverlapp(pred_clusters, go_clusters)
	line = "%s\t%i\t%s\t%i\t%i\t%i\t%i" % (feature_combination, num_training_ppi, "\t".join(map(str, eval_scores)), len(predicted_ppis), len(pred_clusters) , matched_corum_clusters, matched_go_clusters)
	outFH = open("%s.eval.txt" % (outDir, "w"))
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

def main():
	benchmark()
	#cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
