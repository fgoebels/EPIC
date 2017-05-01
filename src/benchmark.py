from __future__ import division
import numpy as np
import copy, os, sys, re
import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils
import matplotlib.pyplot as plt
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

def cut(args):
	fc, scoreF, outF = args
	if fc == "00000000": sys.exit()
	this_scores = get_fs_comb(fc)
	scoreCalc = CS.CalculateCoElutionScores("", "", "","", cutoff=0.5)
	empty_gs = GS.Goldstandard_from_Complexes()
	empty_gs.positive = set([])
	empty_gs.negative = set([])
	scoreCalc.readTable(scoreF, empty_gs)
	print scoreCalc.to_predict
	feature_comb = feature_selector([fs.name for fs in this_scores], scoreCalc, [])
	feature_comb.open()
	outFH = open(outF, "w")
	print >> outFH, "\t".join(feature_comb.scoreCalc.header)
	for i in range(feature_comb.to_predict):
		edge, edge_scores = feature_comb.get_next()
		if edge == "" or edge_scores == []: continue
		print >> outFH, "%s\t%s" % (edge, "\t".join(map(str, edge_scores)))
	outFH.close()
	feature_comb.close()


def merge_MS(args):
	def read_scores(scoreF, cutoff):
		num_prots = CS.lineCount(scoreF)
		scoreFH = open(scoreF)
		header = scoreFH.readline().rstrip()
		header = header.split("\t")
		out = CS.CalculateCoElutionScores("", "", "", 4)
		out.scores = np.zeros((num_prots , len(header[2:])))
		out.header = header
		i = 0
		for line in scoreFH:
			line = line.rstrip()
			if line == "":continue
			line = line.split("\t")
			edge = "\t".join(line[:2])
			this_score = np.array(map(float, line[2:]))
			if len(list(set(np.where(this_score > cutoff)[0]))) > 0:
				out.ppiToIndex[edge] = i
				out.IndexToPpi[i] = edge
				out.scores[i, :] = this_score
				i += 1
		out.scores = out.scores[0:i, :]
		print i
		return out

	ms1_in, ms2_in, mode, outF = args
	ms1cutoff = 0
	ms2cutoff = 0

	if mode == "i":
		ms1cutoff = 0.5
		ms2cutoff = 0.5
	if mode == "u":
		ms1cutoff = 0
		ms2cutoff = 0
	if mode == "l":
		ms1cutoff = 0
		ms2cutoff = 0.5
	if mode == "r":
		ms1cutoff = 0.5
		ms2cutoff = 0

	ms1 = read_scores(ms1_in, ms1cutoff)
	print "Done reading in MS1"
	print ms1.scores.shape

	ms2 = read_scores(ms2_in, ms2cutoff)
	print "Done reading in MS2"
	print ms2.scores.shape


	ms2.merge(ms1, mode)
	print "Done merging MS1 and MS2"
	print ms2.scores.shape


	outFH = open(outF, "w")
	print >> outFH, "\t".join(ms2.header)
	for i in range(ms2.scores.shape[0]):
		if len(list(set(np.where(ms2.scores[i, :] > 0.5)[0]))) > 0:
			print >> outFH, "%s\t%s" % (ms2.IndexToPpi[i], "\t".join(map(str, ms2.scores[i, :])))
	outFH.close()

def exp_comb(args):
	FS, i, j, num_iter, input_dir, num_cores, ref_complexes, scoreF, output_dir = args

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
	if FS == "00000000": sys.exit()
	this_scores = get_fs_comb(FS)
	print this_scores
	no_reference_overlap = False
	fs = False

	#global reference data set, since we don't want to compare with artifical smaller reference data set
	foundprots, _ = utils.load_data(input_dir, [])
	global_gs = Goldstandard_from_cluster_File(ref_complexes, foundprots)

	i,j, num_iter, num_cores = map(int, [i, j, num_iter, num_cores])
	if i == 0 and j == 0: sys.exit()

	out_head = ""
	all_scores = []

	for iter in range(num_iter):
		rnd.seed()
		this_eprofiles = get_eData_comb(input_dir, i, j)
		rnd.seed(1)
		print [f.split(os.sep)[-1] for f in this_eprofiles]
		this_foundprots, _ = utils.load_data(this_eprofiles, [])
		head, scores = run_epic_with_feature_combinations(this_scores, this_eprofiles, num_cores, use_rf, scoreF, output_dir + ".%i_%i.%i" % (i, j, iter),
														  no_reference_overlap, ref_complexes=ref_complexes, fs=fs, globalGS=global_gs)

		out_head = head
		all_scores.append(scores)
		print head
		print scores

	outFH = open(output_dir + ".%i_%i.all.eval.txt" % (i, j), "w")
	print >> outFH, "Num_iex\tNum_beads\t%s" % out_head
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




def EPIC_eval_fs(args):
	in_dir, refF, outF = args
	ref_clusters = GS.Clusters(False)
	ref_clusters.read_file(refF)
	outFH = open(outF, "w")
	i = 0
	allFiles = paths = [os.path.join(in_dir,fn) for fn in next(os.walk(in_dir))[2]]
	for file in allFiles:
		if not file.endswith("clust.txt"): continue
		pred_clusters = GS.Clusters(False)
		pred_clusters.read_file(file)
		_, overlap, _ = pred_clusters.get_matching_complexes(ref_clusters)
		filesplit = file.split(".")[0:4]
		fs_comp = filesplit[0].split(os.sep)[-1]
		scores, head =  utils.clustering_evaluation(ref_clusters, pred_clusters, "Eval", False)
		if i == 0:
			print "FS_code\tCLF\tSE\tFS\tNum_complexes" + "\t".join(np.array(head.split("\t"))[[0,1,6]])
			print >> outFH, "FS_code\tCLF\tSE\tFS\tNum_complexes" + "\t".join(np.array(head.split("\t"))[[0,1,6]])
		print fs_comp + "\t" + "\t".join(filesplit[1:4]) + "\t"+ str(len(pred_clusters.complexes))+"\t" + "\t".join(np.array(scores.split("\t"))[[0, 1, 6]])
		print >> outFH, fs_comp + "\t" + "\t".join(filesplit[1:4]) + "\t"+ str(len(pred_clusters.complexes))+"\t" + "\t".join(np.array(scores.split("\t"))[[0, 1, 6]])
		i += 1
	outFH.close()

def EPIC_eval_fs_DIST(args):

	def getScore(scores):
		out = []
		for cat in [[2], [3], [8]]:
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

	fs_eval_Files, outDir = args

	header, fs, vals, zvals = epic_read_eval(fs_eval_Files)
	fs = np.array(fs)

	def make_hists(x, header, outDir):
		fig = plt.figure()
		plt.rcParams.update({'font.size': 8})
		for i in range(len(header)):
			scores = x[:,i]

			ax = fig.add_subplot(3,4, i+1)
			ax.hist(scores)
			ax.set_title(header[i].replace(" ", "_"), fontsize=6)
			for tick in ax.get_xticklabels():
				tick.set_rotation(45)

		fig.tight_layout()
		fig.subplots_adjust(top=0.88)
		plt.savefig(outDir + ".hist.pdf")
		plt.close()

#	filtering feature selection with low number of predicted clusters

#	sel_vals = set(range(len(vals[:,1])))
#	for k in range(len(header)):
#		this_lb = np.percentile(vals[:, k], 5)
#		this_ub = np.percentile(vals[:, k], 95)
#		sel_vals &= set(np.where(vals[:, k] > this_lb)[0]) & set(np.where(vals[:, k] < this_ub)[0])


	print header[1]
	sel_vals = np.where(vals[:, 1] > 100)[0]
	sel_vals = list(sel_vals)

	vals = vals[sel_vals,]
	zvals = zvals[sel_vals,]
	fs = np.array(fs)[sel_vals,]

	make_hists(vals, header, outDir + ".raw")
	make_hists(zvals, header, outDir + ".zscore")

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


	composit_scores = {}
	scores = {}
	for i in range(len(fs)):
		this_f = fs[i]
		this_vals = getScore(zvals[i,:])
		this_dist = dist(this_vals,max_vals)
		summed_zscores = sum(this_vals)
		composit_scores[this_f] = summed_zscores
		scores[this_f] = this_dist

	scores_sorted = sorted(scores, key=scores.get)

	outFH = open(outDir + ".results.txt", "w")
	print >> outFH, "Artifical optimal vector"
	print >> outFH, "\t" + "\t".join(header)
	print >> outFH, "Scores:\t\t\t" + "\t".join(map(lambda x : "%.2f" % x, this_max_vals))
	print >> outFH, "Optimal FS per category:\t" + "\t".join(this_max_val_fc)

	print >> outFH, "Composit_scorte\tDistance\tFS\t" + "\t".join(header)
	for  f in scores_sorted:
		score = scores[f]
		c_score = composit_scores[f]
		f_scores = "\t".join(map(lambda x : "%.2f" % x, vals[np.where(fs == f)[0], :][0]))
		print >> outFH, "%.2f\t%.2f\t%s\t%s" % (c_score, score, f, f_scores)
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
		self.to_predict = scoreCalc.to_predict
		self.scoreCalc = self.filter_scoreCalc(scoreCalc)


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
		if self.valid_score(scores[self.to_keep_score]):
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
		print "Num val ppis: %i" % len(to_keep_rows)
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[to_keep_rows, :]
		self.update_ppi_map(filtered_scoreCalc, to_keep_rows)
		return filtered_scoreCalc

	def filter_scoreCalc(self, scoreCalc):
		if scoreCalc.scores.shape[0] == 0: # scores are empty no filtering necessary
			filtered_scoreCalc = copy.deepcopy(scoreCalc)
			filtered_scoreCalc.header = np.array(filtered_scoreCalc.header)[self.to_keep_header]
			filtered_scoreCalc.scores = filtered_scoreCalc.scores[:, self.to_keep_score]
			return filtered_scoreCalc
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
			print "recieved empty edge"
			return "", []
		protA, protB = edge.split("\t")

		if (protA not in self.valprots or protB not in self.valprots) and self.valprots != []:
			print "not elution profile for edge %s\t%s" % (protA, protB)
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
	gs = utils.create_goldstandard(taxid, "")
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
		tmp_head, tmp_scores = run_epic_with_feature_combinations(bayes_comb, input_dir, 4, True, scoreF, output_dir ,
										   no_overlap_in_training=False, ref_complexes=ref_compF)
		out_head = tmp_head
		out_scores.append(tmp_scores)

	outFH = open(output_dir + "all.eval.txt", "w")
	print >> outFH, "%s\n%s" % (out_head, "\n".join(out_scores))
	outFH.close()

	print "%s\n%s" % (out_head, "\n".join(out_scores))

def run_epic_with_feature_combinations(feature_combination, input_dir, num_cores, use_rf, scoreF,  output_dir, no_overlap_in_training, ref_complexes="", taxid="", fs = "", globalGS = ""):
	if ref_complexes !="" and taxid!="":
		print "Suplly either taxid or reference complex file"
		sys.exit()
	all_gs = ""

	clf = CS.CLF_Wrapper(num_cores, use_rf, useFeatureSelection=fs)

	foundprots, elution_datas = utils.load_data(input_dir, [])
	if globalGS == "":
		if taxid != "": all_gs = utils.create_goldstandard(taxid, foundprots)
		if ref_complexes != "":
			print "Loading reference from file"
			print ref_complexes
			all_gs = Goldstandard_from_cluster_File(ref_complexes, "")
	else:
		all_gs = globalGS

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

	print len(all_gs.complexes.complexes)
	print len(train.complexes.complexes)
	print len(eval.complexes.complexes)

 	utils.bench_clf(feature_comb, train, eval, clf, "%s.%s" % (output_dir, out_prefix), verbose=True)

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
		if globalGS!="":
			g_s, g_h = utils.clustering_evaluation(globalGS.complexes, pred_clusters, "Global", True)
			out_head += "\t%s" % g_h
			out_scores += "\t%s" % g_s
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
	this_scores = get_fs_comb(feature_combination)
	num_cores = int(num_cores)
	use_rf = use_rf == "True"
	fs = fs =="True"

	no_reference_overlap = False
	head, scores = run_epic_with_feature_combinations(this_scores, input_dir, num_cores, use_rf, scoreF, output_dir,
															  no_reference_overlap, ref_complexes=ref_complexes, fs = fs)
	outFH = open(output_dir + ".eval.wo.txt", "w")
	print >> outFH, "%s\n%s" % (head, scores)
	outFH.close()

def get_fs_comb(comb_string):
	#Create feature combination
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
	this_scores = []
	for i, feature_selection in enumerate(comb_string):
		if feature_selection == "1": this_scores.append(scores[i])
	return this_scores


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

	elif mode == "-merge_ms":
		merge_MS(sys.argv[2:])

	elif mode == "-cut":
		cut(sys.argv[2:])
	elif mode == "-best_fs2":
		EPIC_eval_fs(sys.argv[2:])

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)