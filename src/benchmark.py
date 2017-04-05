from __future__ import division
import numpy as np
import copy, os, sys
import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
from scipy.spatial import distance
from scipy.stats import zscore
import glob
import random as rnd

"""
from scipy.spatial import distance
from scipy.stats import zscore
import itertools

class Goldstandard_from_PPIs():
	def __init__(self, refF, ratio= -1):
		self.goldstandard_positive = set([])
		self.goldstandard_negative = set([])
		allprots = set([])
		refFH = open(refF)
		for line in refFH:
			line = line.rstrip()
			(idA, idB) = line.split("\t")
			allprots.add(idA)
			allprots.add(idB)
			pos_edge = "\t".join(sorted([idA, idB]))
			self.goldstandard_positive.add(pos_edge)
		refFH.close()
		allprots = list(allprots)
		if ratio == -1:
			for i in range(len(allprots)):
				prot_i = allprots[i]
				for j in range(i+1, len(allprots)):
					prot_j = allprots[j]
					edge = "\t".join(sorted([prot_i, prot_j]))
					if edge in self.goldstandard_positive: continue
					self.goldstandard_negative.add(edge)
		else:
			while(len(self.goldstandard_negative)< len(self.goldstandard_positive)*ratio):
				prot_i, prot_j = rnd.sample(allprots,2)
				edge = "\t".join(sorted([prot_i, prot_j]))
				if edge in self.goldstandard_positive: continue
				self.goldstandard_negative.add(edge)
"""

def exp_comb(args):
	i, j, num_iter, input_dir, num_cores, ref_complexes, scoreF, output_dir = args

	def get_eData_comb(data_dir, num_iex, num_beads):
		all_exp =  map(str, glob.glob(data_dir + "*.txt"))
		iex_exp = [f for f in all_exp if (f.split(os.sep)[-1].startswith("all"))]
		beads_exp = [ f for f in all_exp if ( not f.split(os.sep)[-1].startswith("all"))]
		if(i>len(iex_exp)):
			print "i is to large"
			sys.exit()
		if (j > len(beads_exp)):
			print "j is to large"
			sys.exit()

		sel_iex = rnd.sample(iex_exp, num_iex)
		sel_beads = rnd.sample(beads_exp, num_beads)
		return sel_iex + sel_beads


	# EPIC paramters
	use_rf = True
	this_scores = [ CS.MutualInformation(), CS.Bayes(3), CS.Jaccard(), CS.Poisson(5)]
	no_reference_overlap = False
	fs = False

	i,j, num_iter, num_cores = map(int, [i, j, num_iter, num_cores])
	if i == 0 and j == 0: sys.exit()

	out_head = ""
	all_scores = []

	for iter in range(num_iter):
		this_eprofiles = get_eData_comb(input_dir, i, j)
		print this_eprofiles
		head, scores = run_epic_with_feature_combinations(this_scores, this_eprofiles, num_cores, use_rf, scoreF, output_dir + ".%i" % iter,
														  no_reference_overlap, ref_complexes=ref_complexes, fs=fs)

		out_head = head
		all_scores.append(scores)
		print head
		print scores

	outFH = open(output_dir + ".%i_%i.all.eval.txt" % (i, j), "w")
	print >> outFH, "Num_iex,\tNum_beads\t%s" % out_head
	for score in all_scores:
		print >> outFH, "%i\t%i\t%s" % (i,j, score)
	outFH.close()


def EPIC_cor(args):
	print "fuu"
	fs_eval_dir = args[0]
	vals = []
	fs = []
	for fs_eval_F in os.listdir(fs_eval_dir):
		if not fs_eval_F.endswith(".eval.txt"): continue
		fs_eval_F = fs_eval_dir + os.sep + fs_eval_F
		print fs_eval_F
		fs_eval_FH = open(fs_eval_F)
		print fs_eval_FH.readline().split("\t")[1:19]
		for line in fs_eval_FH:
			line = line.rstrip().split("\t")
			if (int(line[2])) < 100: continue
			vals.append(map(float, line[1:19]))
			fs.append(line[0])
		fs_eval_FH.close()
	vals = np.array(vals)
	for row in np.corrcoef(np.transpose(vals)):
		print "\t".join(map("{:.2f}".format, row))





def EPIC_eval_fs_DIST(args):

	def getScore(scores):
		out = []
		for cat in [[0],[1, 2], [3, 4], [6], [7], [9]]:
			out.append(sum(scores[cat])/len(cat))
		return out

	def dist(a,b):
		return distance.euclidean(a,b)

	def epic_read_eval(fs_eval_dir):

		def norm_score(scores, columns):
			for i in columns:
				min_score = min(scores[:, i])
				max_score = max(scores[:, i])
				for j in range(len(scores[:, i])):
					scores[j, i] = (scores[j, i] - min_score) / (max_score - min_score)

		vals = []
		fs = []
		header = ""
		for fs_eval_F in  os.listdir(fs_eval_dir):
			if not fs_eval_F.endswith(".eval.txt"): continue
			fs_eval_F = fs_eval_dir + os.sep + fs_eval_F
			fs_eval_FH = open(fs_eval_F)
			param = "-".join(fs_eval_F.split(os.sep)[-1].split(".")[0:2])
			header = np.array(fs_eval_FH.readline().strip().split("\t"))
			header =  np.append(header[1:3], header[11:])

			print fs_eval_F
			for line in fs_eval_FH:
				line = line.rstrip().split("\t")
				scores  = np.append(line[1:3], line[11:])
				vals.append(map(float, scores))
				fs.append("%s-%s" % (param, line[0]))

			fs_eval_FH.close()

		vals = np.array(vals)
		zvals = zscore(vals)

		return header, fs, vals, zvals

	all_scores = {}

	fs_eval_Files, outfile = args

	header, fs, vals, zvals = epic_read_eval(fs_eval_Files)
	fs = np.array(fs)

#	filtering feature selection with low number of predicted clusters
	sel_vals = list(np.where(vals[:, 1] >= 100)[0])
	vals = vals[sel_vals,]
	zvals = zvals[sel_vals,]
	fs = np.array(fs)[sel_vals,]


	this_max_val_fc = []
	this_max_vals = []
	max_zvals = []
	for i in range(len(header)):
		max_index = np.argmax(np.array(zvals[:, i]))
		this_max_vals.append( vals[max_index, i])
		this_max_val_fc.append(fs[max_index])
		max_zvals.append(zvals[max_index, i])
	max_vals = np.array(getScore(np.array(max_zvals)))

#	print max_zvals
#	print max_vals

	scores = {}
	for i in range(len(fs)):
		this_f = fs[i]
		this_vals = getScore(zvals[i,:])
		this_dist = dist(this_vals,max_vals)
		scores[this_f] = this_dist

	best_dist = min(scores.values())
	scores_sorted = sorted(scores, key=scores.get)

	outFH = open(outfile, "w")
	print >> outFH, "Artifical optimal vector"
	print >> outFH, "\t" + "\t".join(header)
	print >> outFH, "Scores:\t" + "\t".join(map(lambda x : "%.2f" % x, this_max_vals))
	print >> outFH, "Optimal FS per category:\t" + "\t".join(this_max_val_fc)

	print >> outFH, "Distance\tFS\t" + "\t".join(header)
	for  f in scores_sorted:
		score = scores[f]
		f_scores = "\t".join(map(lambda x : "%.2f" % x, vals[np.where(fs == f)[0], :][0]))
		print >> outFH, "%.2f\t%s\t%s" % (score, f, f_scores)
	outFH.close()

def Goldstandard_from_cluster_File(gsF, foundprots = ""):
		clusters = GS.Clusters(need_to_be_mapped=False)
		clusters.read_file(gsF)
		if foundprots != "": clusters.remove_proteins(foundprots)
		gs = GS.Goldstandard_from_Complexes("All")
		gs.complexes = clusters
		gs.make_pos_neg_ppis()
		return gs

class feature_selector:
	def __init__(self, feature_names, scoreCalc, valprots):
		self.valprots = valprots
		self.get_cols(scoreCalc.header, feature_names)
		self.cutoff = scoreCalc.cutoff
		self.scoreCalc = self.filter_scoreCalc(scoreCalc)
		self.to_predict = self.scoreCalc.to_predict

	def set_cutoff(self, cutoff):
		self.cutoff = cutoff

	def get_cols(self, header, feature_names):
		self.to_keep_header = [0, 1]
		self.to_keep_score = []
		for i in range(2, len(header)):
			colname = header[i]
			scorename = colname.split(".")[-1]
			if scorename in feature_names:
				self.to_keep_header.append(i)
				self.to_keep_score.append(i - 2)

	def valid_score(self, scores):
		return len(list(set(np.where(scores > self.cutoff)[0])))>0

	def filter_score(self, scores):
		if self.valid_score(scores):
			return scores[self.to_keep_score]
		else:
			return []

	def update_ppi_map(self, scoreCalc, val_indices):
		prot_indices = [None] * len(scoreCalc.ppiToIndex.keys())
		for p, i in scoreCalc.ppiToIndex.items(): prot_indices[i] = p
		prot_indices = np.array(prot_indices)[val_indices]
		scoreCalc.ppiToIndex = {}
		scoreCalc.IndexToPpi = {}
		for i, p in enumerate(prot_indices):
			scoreCalc.ppiToIndex[p] = i
			scoreCalc.IndexToPpi[i] = p

	def filter_valprots(self, scoreCalc):
		filtered_scoreCalc = copy.deepcopy(scoreCalc)
		to_keep_rows = []
		for p, i in scoreCalc.ppiToIndex.items():
			pA, pB = p.split("\t")
			if pA in self.valprots and pB in self.valprots:
				to_keep_rows.append(i)
		to_keep_rows = np.array(to_keep_rows)
		print "NUm val ppis: %i" % len(to_keep_rows)
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[to_keep_rows, :]
		self.update_ppi_map(filtered_scoreCalc, to_keep_rows)
		return filtered_scoreCalc

	def filter_scoreCalc(self, scoreCalc):
		print scoreCalc.scores.shape
		filtered_scoreCalc = copy.deepcopy(scoreCalc)
		if self.valprots != "": filtered_scoreCalc = self.filter_valprots(filtered_scoreCalc)
		print filtered_scoreCalc.scores.shape
		filtered_scoreCalc.header = np.array(filtered_scoreCalc.header)[self.to_keep_header]
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[:, self.to_keep_score]
		val_rows = np.where(np.apply_along_axis( self.valid_score, 1, filtered_scoreCalc.scores)==True)[0]
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[val_rows, :]
		self.update_ppi_map(filtered_scoreCalc, val_rows)
		print filtered_scoreCalc.scores.shape
		return filtered_scoreCalc

 		return filtered_scoreCalc

	def get_next(self):
		edge, scores = self.scoreCalc.get_next()
		if edge =="":
			return "", []
		protA, protB = edge.split("\t")
		if protA not in self.valprots or protB not in self.valprots:
			return "", []
		return edge, self.filter_score(scores)

	def toSklearnData(self, gs):
		return self.scoreCalc.toSklearnData(gs)

	def open(self):
		self.scoreCalc.open()

	def close(self):
		self.scoreCalc.close()

def write_reference(args):
	input_dir, taxid, output_dir = args
	foundprots, elution_datas = utils.load_data(input_dir, [])
	gs = utils.create_goldstandard(taxid, foundprots)
	out = gs.complexes.to_string()
	outFH = open(output_dir, "w")
	print >> outFH, out
	outFH.close()

def bench_Bayes(args):
	input_dir, scoreF, ref_compF, output_dir = args
	out_head, out_scores = "", []
	combinations = [ [CS.Bayes(1)]
					,[CS.Bayes(2)]
					,[CS.Bayes(3)]]

	for bayes_comb in combinations:
		tmp_head, tmp_scores = run_epic_with_feature_combinations(bayes_comb, input_dir, 4, True, scoreF, output_dir,
										   no_overlap_in_training=False, ref_complexes=ref_compF)
		out_head = tmp_head
		out_scores.append(tmp_scores)

	outFH = open(output_dir + "all.eval.txt", "w")
	print >> outFH, "%s\n%s" % (out_head, "\n".join(out_scores))
	outFH.close()

	print "%s\n%s" % (out_head, "\n".join(out_scores))

def run_epic_with_feature_combinations(feature_combination, input_dir, num_cores, use_rf, scoreF,  output_dir, no_overlap_in_training, ref_complexes="", taxid="", fs = ""):
	if ref_complexes !="" and taxid!="":
		print "Suplly either taxid or reference complex file"
		sys.exit()
	all_gs = ""

	clf = CS.CLF_Wrapper(num_cores, use_rf, useFeatureSelection=fs)

	foundprots, elution_datas = utils.load_data(input_dir, feature_combination)
	if taxid != "": all_gs = utils.create_goldstandard(taxid, foundprots)
	if ref_complexes != "":
		print "Loading reference from file"
		print ref_complexes
		all_gs = Goldstandard_from_cluster_File(ref_complexes, foundprots)

	print len(all_gs.positive)
	print len(all_gs.negative)

	scoreCalc = CS.CalculateCoElutionScores(feature_combination, elution_datas, output_dir + ".scores.txt",
											num_cores=num_cores, cutoff=0.5)
	scoreCalc.readTable(scoreF, all_gs)

	print scoreCalc.scores.shape
	print "Num valid filterd ppis: %i" % len(set(scoreCalc.ppiToIndex.keys()))
	print "Num valid all pos: %i" % len(set(scoreCalc.ppiToIndex.keys()) & set(all_gs.positive))
	print "Num valid all negative: %i" % len(set(scoreCalc.ppiToIndex.keys()) & set(all_gs.negative))
	all_gs.make_pos_neg_ppis(val_ppis=set(scoreCalc.ppiToIndex.keys()))
	feature_comb = feature_selector([fs.name for fs in feature_combination], scoreCalc, foundprots)
	print feature_comb.scoreCalc.scores.shape

	print "Num valid filterd ppis: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()))
	print "Num valid filtered pos: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()) & set(all_gs.positive))
	print "Num valid filtered negative: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()) & set(all_gs.negative))

	out_prefix = "_".join([fs.name for fs in feature_combination])
	train, eval = all_gs.split_into_holdout_training(set(feature_comb.scoreCalc.ppiToIndex.keys()), no_overlapp=no_overlap_in_training)

# 	utils.bench_clf(feature_comb, train, eval, clf, "%s.%s" % (output_dir, out_prefix), verbose=True)

	print "Num valid ppis in training pos: %i" % len(train.positive)
	print "Num valid ppis in training neg: %i" % len(train.negative)

	print "Num valid ppis in eval pos: %i" % len(eval.positive)
	print "Num valid ppis in eval neg: %i" % len(eval.negative)

	network = utils.make_predictions(feature_comb, "exp", clf, train, "", verbose=True)

	outFH = open("%s.%s.pred.txt" % (output_dir, out_prefix), "w")
	print >> outFH, "\n".join(network)
	outFH.close()

	num_ppis = len(network)
	if num_ppis != 0:

		# Predicting clusters
		utils.predict_clusters("%s.%s.pred.txt" % (output_dir, out_prefix),
							   "%s.%s.clust.txt" % (output_dir, out_prefix))

		# Evaluating predicted clusters
		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file("%s.%s.clust.txt" % (output_dir, out_prefix))
		num_cluster = CS.lineCount("%s.%s.clust.txt" % (output_dir, out_prefix))

		t_s, t_h = utils.clustering_evaluation(train.complexes, pred_clusters, "Train", True)
		e_s, e_h = utils.clustering_evaluation(eval.complexes, pred_clusters, "Eval", True)
		out_head = "FS\tNum PPIs\tNum Clusters\t%s\t%s" % (t_h, e_h)
		out_scores = "%s\t%i\t%i\t%s\t%s" % (out_prefix, num_ppis, num_cluster, t_s, e_s)
		return out_head, out_scores
	else:
		head_t = "\t".join(["%s %s" % ("Train", h) for h in
						  ["mmr", "overlapp", "simcoe", "mean_simcoe_overlap", "sensetivity", "ppv", "accuracy",
						   "sep"]])
		head_e = "\t".join(["%s %s" % ("Eval", h) for h in
						  ["mmr", "overlapp", "simcoe", "mean_simcoe_overlap", "sensetivity", "ppv", "accuracy",
						   "sep"]])

		return "FS\tNum PPIs\tNum Clusters\t%s\t%s" % (head_t, head_e), "\t".join(map(str, [0]*18))

def calc_feature_combination(args):
	feature_combination, input_dir, use_rf, fs, num_cores, scoreF, ref_complexes, output_dir = args
	#Create feature combination
	if feature_combination == "00000000": sys.exit()
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(feature_combination):
		if feature_selection == "1": this_scores.append(scores[i])
	num_cores = int(num_cores)
	use_rf = use_rf == "True"
	fs = fs =="True"

	no_reference_overlap = False
	head, scores = run_epic_with_feature_combinations(this_scores, input_dir, num_cores, use_rf, scoreF, output_dir,
															  no_reference_overlap, ref_complexes=ref_complexes, fs = fs)
	outFH = open(output_dir + ".eval.wo.txt", "w")
	print >> outFH, "%s\n%s" % (head, scores)
	outFH.close()

def main():
	mode = sys.argv[1]

	if mode == "-fs":
		calc_feature_combination(sys.argv[2:])

	elif mode == "-make_ref":
		write_reference(sys.argv[2:])

	elif mode == "-bench_bayes":
		bench_Bayes(sys.argv[2:])

	elif mode == "-cor_eval":
		EPIC_cor(sys.argv[2:])

	elif mode == "-best_fs":
		EPIC_eval_fs_DIST(sys.argv[2:])

	elif mode == "-exp_comb":
		exp_comb(sys.argv[2:])


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)