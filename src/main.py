from __future__ import division

import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import sys

def main():
	feature_combination, input_dir, use_rf, num_cores, mode, anno_source, target_taxid, output_dir = sys.argv[1:]

	#Create feature combination
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])

	print "\t".join([fs.name for fs in this_scores])

	# Initialize CLF
	use_rf = use_rf == "True"
	num_cores = int(num_cores)
	clf = CS.CLF_Wrapper(num_cores, use_rf)

	# Load elution data
	foundprots, elution_datas = utils.load_data(input_dir, this_scores)

	# Generate reference data set
	all_gs = utils.create_goldstandard(target_taxid, foundprots)

#	print len(all_gs.positive)
#	print len(all_gs.negative)


	scoreCalc = CS.CalculateCoElutionScores(this_scores, elution_datas, output_dir + ".scores.txt", num_cores=num_cores, cutoff= 0.5)
#	scoreCalc.calculate_coelutionDatas(all_gs)
	scoreCalc.readTable(output_dir + ".scores.txt", all_gs)
	print len(set(scoreCalc.ppiToIndex.keys()))
	train, eval = all_gs.split_into_holdout_training(set(scoreCalc.ppiToIndex.keys()))

	# Evaluate classifier
	utils.bench_clf(scoreCalc, train, eval, clf, output_dir, verbose=True)

	functionalData = ""
	if mode != "exp":
		functionalData = utils.get_FA_data(anno_source)


	# Predict protein interaction
	network =  utils.make_predictions(scoreCalc, mode, clf, train, functionalData)
	outFH = open("%s.%s.pred.txt" % (output_dir, mode), "w")
	print >> outFH, "\n".join(network)
	outFH.close()

	# Predicting clusters
	utils.predict_clusters("%s.%s.pred.txt" % (output_dir, mode), "%s.%s.clust.txt" % (output_dir, mode))

	# Evaluating predicted clusters
	pred_clusters = GS.Clusters(False)
	pred_clusters.read_file("%s.%s.clust.txt" % (output_dir, mode))
#	utils.clustering_evaluation(train.complexes, pred_clusters, "Train", True)
	utils.clustering_evaluation(eval.complexes, pred_clusters, "", True)


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)
