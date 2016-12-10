from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import numpy as np
import sys
import os

from itertools import chain, combinations


def calculate_allscores():
	number_of_cores, elutionFiles, outDir = sys.argv[1:]
	number_of_cores = int(number_of_cores)

	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	_, elution_data = CS.load_data(elutionFiles, scores)

	scoreCalc = CS.CalculateCoElutionScores()
	scoreCalc.calculate_coelutionDatas(elution_data, scores, outDir, number_of_cores)

def benchmark():
	feature_combination, use_random_forest, number_of_cores, elutionFiles, all_scoreF, outDir = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
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
#	scoreCalc.readTable(all_scoreF)
	scoreCalc.calculate_coelutionDatas(elution_datas, this_scores, outDir, number_of_cores)

	scoreCalc.addLabels(all_p, all_n)
	CS.predictInteractions(scoreCalc, outDir , use_random_forest, number_of_cores)
	predF = "%s.pred.txt" % (outDir)
	clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	os.system(clustering_CMD)

	scoreCalc.addLabels(training_p, training_n)
#	CS.bench_scores(scoreCalc, outDir, number_of_cores, useForest=use_random_forest)
	CS.clustering_evaluation(evals, scoreCalc, outDir, ",".join([score.name for score in this_scores]), number_of_cores, use_random_forest)

def main():
	calculate_allscores()
	#benchmark()
	#cluster_overlapp()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
