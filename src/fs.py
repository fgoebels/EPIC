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
				out.append("%s\t%f" % (edge, pred_prob))
	All_score_FH.close()
	return out

def benchmark2():
	feature_combination, use_random_forest, number_of_cores, elutionFiles, all_scoreF, outDir = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(50), CS.Pearson(), CS.Apex()]
	this_scores = []
	number_of_cores = int(number_of_cores)
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])
	print this_scores

	foundprots, elution_datas = CS.load_data(elutionFiles, this_scores)
	evals = CS.create_goldstandard("6239", foundprots)
	training_p, training_n, all_p, all_n, go_complexes, corum_complexes = evals
	scoreCalc = CS.CalculateCoElutionScores()
	scoreCalc.addLabels(all_p, all_n)
	scoreCalc.readTable(all_scoreF)
#	scoreCalc.calculate_coelutionDatas(elution_datas, this_scores, outDir, number_of_cores)

	scoreCalc.addLabels(all_p, all_n)
	CS.predictInteractions(scoreCalc, outDir , use_random_forest, number_of_cores)
	predF = "%s.pred.txt" % (outDir)
	clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	os.system(clustering_CMD)

	scoreCalc.addLabels(training_p, training_n)
	#CS.bench_scores(scoreCalc, outDir, number_of_cores, useForest=use_random_forest)
	CS.clustering_evaluation(evals, scoreCalc, outDir, ",".join([score.name for score in this_scores]), number_of_cores, use_random_forest)

def benchmark():
	feature_combination, use_random_forest, number_of_cores, elutionFiles, refF, train_scoreF, all_scoreF, go_complexesF, corum_complexesF, outDir = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(50), CS.Pearson(), CS.Apex()]
	use_random_forest = use_random_forest == True
	number_of_cores = int(number_of_cores)
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])

	print this_scores

	elutionFH = open(elutionFiles)
	elutionDatas = []
	elutionProts = set([])
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = CS.ElutionData(elutionFile)
		elutionDatas.append(elutionData)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())

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

	go_cluster = GS.Goldstandard_from_cluster_File(go_complexesF, found_prots=elutionProts)
	go_cluster = go_cluster.clusters
	print len(go_cluster.complexes)

	corum_cluster = GS.Goldstandard_from_cluster_File(corum_complexesF, found_prots=elutionProts)
	corum_cluster = corum_cluster.clusters
	print len(corum_cluster.complexes)

	scorecalc_train.addLabels(positive, negative)
	scorecalc_train.rebalance()

	_, data, targets = scorecalc_train.toSklearnData(get_preds=False)
	num_training_ppi =  data.shape[0]
	data = np.array(data)
	clf = CS.CLF_Wrapper(data, targets, num_cores=number_of_cores, forest=use_random_forest, folds=2, useFeatureSelection=False)
	eval_scores = clf.getValScores()
	print "Done doing ML"
#	predicted_ppis = predictInteractions(all_scoreF, scorecalc_train, this_scores, useForest = use_random_forest, num_cores = number_of_cores)
	predF = "%s.pred.txt" % (outDir)
#	predFH = open(predF, "w")
#	predFH.write("\n".join(predicted_ppis))
#	predFH.close()
	predicted_ppis = CS.lineCount(predF)
	print "Done with prediction"
	clustering_CMD = "java -jar /Users/florian/workspace/EPIC/src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	os.system(clustering_CMD)
	pred_clusters = GS.Clusters(need_to_be_mapped=False)
	pred_clusters.read_file("%s.clust.txt" % (outDir))
	pred_clusters.filter_complexes()
	pred_clusters = pred_clusters
	corum_scores = "\t".join(map(str, pred_clusters.clus_eval(corum_cluster)))
	go_scores = "\t".join(map(str, pred_clusters.clus_eval(go_cluster)))
	line = "%s\t%i\t%s\t%i\t%i\t%s\t%s" % ("\t".join(list(feature_combination)), num_training_ppi, "\t".join(map(str, eval_scores)), predicted_ppis, len(pred_clusters.complexes) , corum_scores, go_scores)
	outFH = open("%s.eval.txt" % (outDir), "w")
	print line
	print >> outFH, line
	outFH.close()

def removeClusters(compA, compB):
	out = {}
	protsB = set([])
	for comp in compB:
		protsB |= compB[comp]
	for comp in compA:
		if len(compA[comp] & protsB) >= 2:
			out[comp] = compA[comp]
	return out

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
	benchmark2()
	#cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
