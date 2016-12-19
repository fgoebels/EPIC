from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import numpy as np
import sys
import os

from itertools import chain, combinations

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

def readTable(feature_scores, scoreF, gs, score_cutoff = 0.5 ):
	out = CS.CalculateCoElutionScores()
	with open(scoreF) as scoreFH:
		header = scoreFH.readline().rstrip().split("\t")
		new_header, scores_to_keep = get_cols(header, feature_scores)
		out.header = new_header
		out.scores = np.zeros((len(gs), len(new_header) - 2))
		i = 0
		out.ppiToIndex = {}
		for line in scoreFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			edge = "\t".join(sorted(linesplit[0:2]))
			if gs !="" and edge not in gs: continue
			edge_scores = np.nan_to_num(np.array(map(float, np.array(linesplit[2:])[scores_to_keep],)))
			if len(list(set(np.where(edge_scores > score_cutoff)[0])))>0:
				out.scores[i, :] = edge_scores
				out.IndexToPpi[i] = edge
				out.ppiToIndex[edge] = i
				i += 1
	print i
	out.scores = out.scores[:i,:] # remove empty zero rows at the end of the matrix
	print out.scores[10,:]
	print out.scores.shape
	scoreFH.close()
	return out, scores_to_keep

def calc_chunkscors():
	number_of_cores, elutionFiles, chunkF  = sys.argv[1:]
	number_of_cores = int(number_of_cores)
	ppis = set([])
	chunkFH = open(chunkF)
	for line in chunkFH:
		line = line.rstrip()
		ppis.add(line)
	chunkFH.close()

	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(),
		  CS.Apex()]
	_, elution_data = CS.load_data(elutionFiles, scores)


	scoreCalc = CS.CalculateCoElutionScores()
	scoreCalc.calculate_coelutionDatas(elution_data, scores, "", number_of_cores, toPred=ppis)
	outFH = open(chunkF+".scores.txt","w")
	scoreCalc.toTable(fh=outFH, labels=False)
	outFH.close()

def calculate_allscores():
	number_of_cores, elutionFiles, chunkF, outDir = sys.argv[1:]
	number_of_cores = int(number_of_cores)
	ppis = set([])
	chunkF

	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	_, elution_data = CS.load_data(elutionFiles, scores)

	scoreCalc = CS.CalculateCoElutionScores()
	scoreCalc.calculate_coelutionDatas(elution_data, scores, outDir, number_of_cores)



def predictInteractions(All_score_F, outDir, train_scoreCalc, scores, useForest = True, num_cores =4, score_cutoff=0.5):
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
	outFH = open(outDir + ".pred.txt", "w")
	outFH.write("\n".join(out))
	outFH.close()

def make_ref_data():
	target_taxid, elution_profiles_dir, out_dir = sys.argv[1:]
	foundprots, elution_datas = CS.load_data(elution_profiles_dir, [])

	training_p, training_n, _, _, go_complexes, corum_complexes = CS.create_goldstandard(target_taxid, foundprots)
	corum_complexes.write_cuslter_file(out_dir + "/corum.comp.txt")
	go_complexes.write_cuslter_file(out_dir + "/go.comp.txt")
	gsFH = open(out_dir + "/gs.txt", "w")
	for ppi in training_p:
		print >> gsFH, "%s\tpositive" % ppi
	for ppi in training_n:
		print >> gsFH, "%s\tnegative" % ppi
	gsFH.close()

	scorecalc = CS.CalculateCoElutionScores()
	topred = scorecalc.getAllPairs_coelutionDatas(elution_datas)
	predFH = open(out_dir + "/topred.txt", "w")
	predFH.write("\n".join(topred))
	predFH.close()


def benchmark():
	feature_combination, use_random_forest, number_of_cores, gsF, corumF, goF, evalF, all_scoreF, outDir = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(50), CS.Pearson(), CS.Apex()]
	this_scores = []
	number_of_cores = int(number_of_cores)
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])
	print this_scores

#	foundprots, elution_datas = CS.load_data(data_dir, this_scores)

	lb=2
	go_complexes = GS.Clusters(False, lb=lb)
	go_complexes.read_file(goF)
	go_complexes.filter_complexes()
	go_complexes.merge_complexes()
	go_complexes.filter_complexes()
#	go_complexes.remove_proteins(foundprots)

	corum_complexes = GS.Clusters(False, lb=lb)
	corum_complexes.read_file(corumF)
	print len(corum_complexes.complexes)
	corum_complexes.filter_complexes()
	corum_complexes.merge_complexes()
	corum_complexes.filter_complexes()
#	corum_complexes.remove_proteins(foundprots)
	print len(corum_complexes.complexes)

	gs = GS.Goldstandard_from_reference_File(gsF)
	training_p, training_n = gs.goldstandard_positive, gs.goldstandard_negative
	print "Done loading gs"

#	Use precalcutaed files since cinet as no internet acces on computational nodes
#	evals = CS.create_goldstandard("6239", foundprots)
#	training_p, training_n, all_p, all_n, go_complexes, corum_complexes = evals

	print all_scoreF
	scoreCalc, scores_to_keep = readTable(this_scores, evalF, gs = training_p | training_n)
	scoreCalc.addLabels(training_p, training_n)
	scoreCalc.rebalance(ratio=5)
	print len(scoreCalc.positive)
	print len(scoreCalc.negative)
	ids_train, data_train, targets_train = scoreCalc.toSklearnData(get_preds=False)
	print data_train.shape
	#predictInteractions(all_scoreF, outDir, scoreCalc, this_scores, use_random_forest, number_of_cores, 0.5)
	CS.predictInteractions(scoreCalc, outDir , use_random_forest, number_of_cores, scoreF=all_scoreF, verbose=True, fs = scores_to_keep)
	predF = "%s.pred.txt" % (outDir)
	clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	os.system(clustering_CMD)

#	scoreCalc.addLabels(training_p, training_n)
#	CS.bench_scores(scoreCalc, outDir, number_of_cores, useForest=use_random_forest)
	CS.clustering_evaluation([training_p, training_n, go_complexes, corum_complexes], scoreCalc, outDir, ",".join([score.name for score in this_scores]), number_of_cores, use_random_forest)

def main():
	#make_ref_data()
	calc_chunkscors()
	#calculate_allscores()
	#benchmark()
	#cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
