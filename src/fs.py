from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import numpy as np
import sys
import os
import glob

from itertools import chain, combinations


def get_cols(header, feature_scores, elution_files = ""):
	scoreNames = set([])
	for st in feature_scores: scoreNames.add(st.name)
	to_keep_header = [0, 1]
	to_keep_score = []
	for i in range(2, len(header)):
		colname = header[i]
		scorename = colname.split(".")[-1]
		filename = colname.rsplit(".", 1)[0]

		if (scorename in scoreNames and elution_files =="") or (scorename in scoreNames and filename in elution_files):
			to_keep_header.append(i)
			to_keep_score.append(i - 2)
	new_header = list(np.array(header)[to_keep_header])
	return new_header, to_keep_score

def readTable(feature_scores, scoreF, gs, score_cutoff = 0.5, outF = "", elution_files = "" ):
	out = CS.CalculateCoElutionScores()
	if outF !="":  outFH = open(outF, "w")
	with open(scoreF) as scoreFH:
		header = scoreFH.readline().rstrip().split("\t")
		new_header, scores_to_keep = get_cols(header, feature_scores, elution_files)
		out.header = new_header
		if outF !="":
			print >> outFH, "\t".join(new_header)
			outFH.flush()
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
				if outF!="":
					print >> outFH, "%s\t%s" % (edge, "\t".join(map(str, edge_scores)))
					continue
				out.scores[i, :] = edge_scores
				out.IndexToPpi[i] = edge
				out.ppiToIndex[edge] = i
				i += 1
	out.scores = out.scores[:i,:] # remove empty zero rows at the end of the matrix
	scoreFH.close()
	if outF != "": outFH.close()
	return out, scores_to_keep

def cut_table():
	feature_combination, all_scoreF, outF = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])
	print this_scores
	readTable(this_scores, all_scoreF,  outF = outF, gs="",)


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

def make_ref_data():
	target_taxid, elution_profiles_dir, out_dir = sys.argv[1:]
	foundprots, elution_datas = CS.load_data(elution_profiles_dir, [])

	training_p, training_n, all_p, all_n, holdout_complexes, training_complexes, all_complexes = CS.create_goldstandard(target_taxid, foundprots)
	holdout_complexes.write_cuslter_file(out_dir + "/holdout.comp.txt")
	all_complexes.write_cuslter_file(out_dir + "/all.comp.txt")
	training_complexes.write_cuslter_file(out_dir + "/training.comp.txt")
	gsFH = open(out_dir + "/train_gs.txt", "w")
	for ppi in training_p:
		print >> gsFH, "%s\tpositive" % ppi
	for ppi in training_n:
		print >> gsFH, "%s\tnegative" % ppi
	gsFH.close()

	gsFH = open(out_dir + "/all_gs.txt", "w")
	for ppi in all_p:
		print >> gsFH, "%s\tpositive" % ppi
	for ppi in all_n:
		print >> gsFH, "%s\tnegative" % ppi
	gsFH.close()

#	scorecalc = CS.CalculateCoElutionScores()
#	topred = scorecalc.getAllPairs_coelutionDatas(elution_datas)
#	predFH = open(out_dir + "/topred.txt", "w")
#	predFH.write("\n".join(topred))
#	predFH.close()


def benchmark(FA, mode):
	feature_combination, use_random_forest, number_of_cores, gs_train_F, gs_all_F, trainingF, holdoutF, allF, eval_scoreF, all_scoreF, outDir = sys.argv[1:]
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	use_random_forest = use_random_forest == "True"
	number_of_cores = int(number_of_cores)
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])
	print ("score is coming")
	print this_scores

#	foundprots, elution_datas = CS.load_data(data_dir, this_scores)

	holdout_comp = GS.Clusters(False)
	holdout_comp.read_file(holdoutF)

	train_comp = GS.Clusters(False)
	train_comp.read_file(trainingF)

	all_comp = GS.Clusters(False)
	all_comp.read_file(allF)

	train_gs = GS.Goldstandard_from_reference_File(gs_train_F)
	training_p, training_n = train_gs.goldstandard_positive, train_gs.goldstandard_negative

	all_gs = GS.Goldstandard_from_reference_File(gs_all_F)
	all_p, all_n = all_gs.goldstandard_positive, all_gs.goldstandard_negative

	header = ""
	line = ""

	#gm = CS.Genemania("6239")
	#gm_scoreC = gm.scoreCalc


	#make_predictions(this_scores, eval_scoreF, all_scoreF, training_p,  training_n, number_of_cores, use_random_forest, outDir + ".train", gm.scoreCalc, mode)
	#sys.exit() # itnegrate Fa without bugs then continue


	#tmp_line, tmp_head = CS.clustering_evaluation(train_comp, "Train", outDir + ".train")
	#train_num_ppis = CS.lineCount(outDir + ".train.pred.txt")
	#train_num_comp = CS.lineCount(outDir + ".train.clust.txt")
	#line += "%s\t%i\t%i" % (feature_combination, train_num_ppis, train_num_comp)
	#header += "Train num pred PPIs\tTrain num pred clust"
	#line += "\t" + tmp_line
	#header += "\t" + tmp_head
	#tmp_line, tmp_head = CS.clustering_evaluation(holdout_comp, "Holdout", outDir + ".train")
	#line += "\t" + tmp_line
	#header += "\t" + tmp_head

	# if FA works make initial eval and see it works
	# print header + "\n" + line
	# sys.exit()



	clusterFileoutDir = make_predictions(this_scores, eval_scoreF, all_scoreF, all_p,  all_n, number_of_cores, use_random_forest, outDir, FA, mode)
	tmp_line, tmp_head = CS.clustering_evaluation(all_comp, "All", clusterFileoutDir[1])
	all_num_ppis = CS.lineCount(clusterFileoutDir[0])
	all_num_comp = CS.lineCount(clusterFileoutDir[1])
	line += "\t%i\t%i" % (all_num_ppis, all_num_comp)
	header += "\tAll num pred PPIs\tAll num pred clust"
	line += "\t" + tmp_line
	header += "\t" + tmp_head



	outFH = open(outDir + ".eval.txt", "w")
	print >> outFH, line
	outFH.close()
	line = line.split("\t")
	header = header.split("\t")
	for i in range(len(line)):
		print "%s\t%s" % (header[i], line[i])

def make_predictions(fc, train_scoreF, all_scoreF, pos, neg, num_cores, use_rf, outDir, fun_anno, mode):
	scoreCalc, scores_to_keep = readTable(fc, train_scoreF, gs=(pos | neg))

	scoreCalc.addLabels(pos, neg)
	scoreCalc.rebalance(ratio=5)
	fun_anno.addLabels(pos,neg)

	#instead of predicting only experimental predict: Exp, FA, and EXP+FA
	#scoreCalc.merge_singe_ScoreCalc(FA_scorecalc)
	#gene name to file
	#genemane.scoreclacl.toString()

	if mode == "exp":
	#predicts using experiment only

		CS.predictInteractions(scoreCalc, outDir + ".exp" , use_rf, num_cores, scoreF=all_scoreF, verbose=True, fs = scores_to_keep)
		outDir = outDir + ".exp"


	#predicts using fun_anno only
	# either switch file rad in memory FA data
	# or write FA data as file and use it as input
	if mode == "fa":

		CS.predictInteractions(fun_anno, outDir + ".fa", use_rf, num_cores, scoreF=all_scoreF, fun_anno_toadd=fun_anno, verbose=True, fs=scores_to_keep, mode="fa", score_cutoff=0)
		outDir = outDir + ".fa"



	#predict using both functional annotation and exp
	# merge scorecalc obj
	# and add FA scores to prediction scores
	if mode == "merge":
		scoreCalc.merge_singe_ScoreCalc(fun_anno,"l")
		outDir = outDir + ".merge"
		CS.predictInteractions(scoreCalc, outDir , use_rf, num_cores, scoreF=all_scoreF, fun_anno_toadd=fun_anno, verbose=True, fs = scores_to_keep, mode="merge")

	# collect the three predicted networks and do merging operation: (EXP union EXP_FA) - (FA - (EXP_FA)) # if the networks are sets you can write (exp | exp_fa) - (fa - exp_fa)
	# be careful to not lose scoring for machine learning method for example EXP predicts A\tB\tS1 and EXP_FA predicts A\tB\tS2 take score S1


	#predict protein complexes using PPIs excluding PPIs only from functional evidence (merged + exp - FA)
	if mode == "final":
		#predict PPIs using only exp data:
		CS.predictInteractions(scoreCalc, outDir + ".exp", use_rf, num_cores, scoreF=all_scoreF, verbose=True,
							   fs=scores_to_keep)
		outDirExp = outDir + ".exp.pred.txt"

		#predict PPIs using Functional evidence only
		CS.predictInteractions(fun_anno, outDir + ".fa", use_rf, num_cores, scoreF=all_scoreF, fun_anno_toadd=fun_anno,
						   verbose=True, fs=scores_to_keep, mode="fa", score_cutoff=0)
		outDirFA = outDir + ".fa.pred.txt"

		#predict PPIs using merged data: functional evidence plus exp
		scoreCalc.merge_singe_ScoreCalc(fun_anno, "l")
		CS.predictInteractions(scoreCalc, outDir + ".merge", use_rf, num_cores, scoreF=all_scoreF, fun_anno_toadd=fun_anno,
						   verbose=True, fs=scores_to_keep, mode="merge")
		outDirMerged = outDir + ".merge.pred.txt"

		finalPPIsDict = {}

		functionalEvidencePPIsDict = {}

		#read exp based PPIs into finalPPIsDict
		with open(outDirExp) as fp:
			for line in fp:
				proteinA, proteinB, score = line.split()
				edge = "\t".join(sorted([proteinA, proteinB]))
				if edge not in finalPPIsDict:
					finalPPIsDict[edge] = score

		#read functional evidence based PPIs into functionalEvidencePPIsDict
		with open(outDirFA) as fp:
			for line in fp:
				proteinA, proteinB, score = line.split()
				edge = "\t".join(sorted([proteinA, proteinB]))
				if edge not in finalPPIsDict:
					functionalEvidencePPIsDict[edge] = score

		#read merged PPIs and if the PPIs only comes from functional evidence, then delete it.
		with open(outDirMerged) as fp:
			for line in fp:
				proteinA, proteinB, score = line.split()
				edge = "\t".join(sorted([proteinA, proteinB]))
				if edge not in functionalEvidencePPIsDict:
					finalPPIsDict[edge] = score

		outDir = outDir + ".final"
		outFH = open(outDir + ".pred.txt", "w")

		for key in finalPPIsDict:
			outFH.write(key)
			outFH.write("\t")
			outFH.write(finalPPIsDict[key])
			outFH.write("\n")
			#eachLine = "\t".join(key, finalPPIsDict[key])
			#outFH.write("\n".join(eachLine))

		outFH.close()

	#predict protein complexes from the PPIs file using CLusterOne algorithm...
	predF = "%s.pred.txt" % (outDir)
	clustering_CMD = "java -jar /Users/lucasminghu/Desktop/EPIC_08_02_2017/EPIC/src/ClusterOne/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	print (clustering_CMD)
	os.system(clustering_CMD)
	print ("debugging here!")
	print (outDir + ".pred.txt")

	return (outDir + ".clust.txt", outDir + ".pred.txt")



def main():
#	cut_table()
#	make_ref_data()
	gm = CS.Genemania("6239")
	#print gm
	sc = gm.getScoreCalc()
	#print sc
	#print sc.IndexToPpi
	#sys.exit()
	benchmark(sc, "final")
#	cluster_overlapp()
#	calc_chunkscors()
#	calculate_allscores()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)
