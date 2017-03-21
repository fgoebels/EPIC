from __future__ import division
import numpy as np
import copy, os, sys
import CalculateCoElutionScores as CS
import GoldStandard as GS
import utils as utils


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


def epic_read_eval(fs_eval_Fs):

	def norm_score(scores, columns):
		for i in columns:
			min_score = min(scores[:,i])
			max_score = max(scores[:,i])
			for j in range(len(scores[:,i])):
				scores[j, i]  = (scores[j,i]-min_score)/(max_score-min_score)

	vals= []
	fs = []
	header = ""
	for fs_eval_F in fs_eval_Fs:
		fs_eval_FH = open(fs_eval_F)
		param = "-".join(fs_eval_F.split(os.sep)[-1].split(".")[0:2])
		header = fs_eval_FH.readline().strip().split("\t")[19:]

		print fs_eval_F
		for line in fs_eval_FH:
			line = line.rstrip().split("\t")
			vals.append(map(float, line[19:]))
			fs.append("%s-%s" % (param, line[0]))

		fs_eval_FH.close()


	vals = np.array(vals)
	zvals = zscore(vals)

	return header, fs, vals, zvals

# Class for selecting sub ser of features of a scoreCalc Object
class feature_selector():
	def init(self):
		return None

def EPIC_eval_fs_DIST():

	def getScore(scores):
		out = []
		for cat in [[0], [1, 2, 9], [3,4,5,6], [7]]:
			out.append(sum(scores[cat])/len(cat))
		return out
	def dist(a,b):
		return distance.euclidean(a,b)

	all_scores = {}

	fs_eval_Files = sys.argv[1:]

	header, fs, vals, zvals = epic_read_eval(fs_eval_Files)
	fs = np.array(fs)

#	filtering feature selection with low number of predicted clusters
#	sel_vals = list(np.where(vals[:, 1] >= 100)[0])
#	vals = vals[sel_vals,]
#	zvals = zvals[sel_vals,]
#	fs = np.array(fs)[sel_vals,]


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
	for  f in scores:
		if scores[f]==best_dist:
			print f
	 		for i in range(len(header)):
				print "%s\t%s\t%s\t%s" % (header[i], str(this_max_vals[i]), this_max_val_fc[i], str(vals[np.where(fs == f)[0], i][0]))


def EPIC_cor():
	fs_eval_dir = sys.argv[1]
	vals = []
	fs = []
	for fs_eval_F in os.listdir(fs_eval_dir):
		if not fs_eval_F.endswith(".eval.txt"): continue
		fs_eval_F = fs_eval_dir + os.sep + fs_eval_F
		print fs_eval_F
		fs_eval_FH = open(fs_eval_F)
		fs_eval_FH.readline()
		for line in fs_eval_FH:
			line = line.rstrip().split("\t")
			if (int(line[2])) < 100: continue
			vals.append(map(float, line[1:19]))
	#		if (int(line[20])) < 100: continue
	#		vals.append(map(float, line[19:]))
			fs.append(line[0])
		fs_eval_FH.close()
	vals = np.array(vals)
	for row in np.corrcoef(np.transpose(vals)):
		print "\t".join(map("{:.2f}".format, row))


def exp_comb():
	i, j, expdir, gsF, corumF, goF, ref_scores_F, all_scores_F, outDir = sys.argv[1:]
	# EPIC paramters
	num_cores = 4
	rf = True
#	this_fs = [CS.Poisson(5), CS.MutualInformation()]
	this_fs = [ CS.Wcc(), CS.Poisson(5)]

	i = int(i)
	j = int(j)

	all_exp =  map(str, glob.glob(expdir + "*.txt"))
	iex_exp = [f for f in all_exp if (f.split(os.sep)[-1].startswith("all"))]
	beads_exp = [ f for f in all_exp if ( not f.split(os.sep)[-1].startswith("all"))]

	if(i>len(iex_exp)):
		print "i is to large"
		sys.exit()
	if (j > len(beads_exp)):
		print "j is to large"
		sys.exit()

	if i == 0 and j == 0: sys.exit()
	out_lines = []
	for k in range(20):
		this_iex   = re.sample(iex_exp,i)
		this_beads = re.sample(beads_exp,j)
		print this_beads
		elution_profiels = this_iex + this_beads
		elutionDatas = []
		foundprots = set([])
		fnames = set([])
		for elutionFile in elution_profiels:
			fnames.add(elutionFile.split(os.sep)[-1])
			elutionData = CS.ElutionData(elutionFile)
			elutionDatas.append(elutionData)
			foundprots = foundprots | set(elutionData.prot2Index.keys())
			for score in this_fs:
				score.init(elutionData)

		go_complexes = GS.Clusters(False)
		go_complexes.read_file(goF)
		go_complexes.remove_proteins(foundprots)
		go_complexes.filter_complexes()

		corum_complexes = GS.Clusters(False)
		corum_complexes.read_file(corumF)
		corum_complexes.remove_proteins(foundprots)
		corum_complexes.filter_complexes()

		gs = GS.Goldstandard_from_reference_File(gsF)
		training_p, training_n = gs.goldstandard_positive, gs.goldstandard_negative
		print "Done loading gs"
#				evals = CS.create_goldstandard(taxid, foundprots)
#				training_p, training_n, all_p, all_n, go_complexes, corum_complexes = evals
#				scoreCalc = CS.CalculateCoElutionScores()
#				scoreCalc.addLabels(all_p, all_n)
#				scoreCalc.calculate_coelutionDatas(elutionDatas, this_scores, output_dir, num_cores)
		scoreCalc, scores_to_keep = fs.readTable(this_fs, ref_scores_F, gs=training_p | training_n, elution_files=fnames)

		scoreCalc.addLabels(training_p, training_n)
		scoreCalc.rebalance(ratio=5)
		print len(scoreCalc.positive)
		print len(scoreCalc.negative)
		print scoreCalc.scores.shape
		ids_train, data_train, targets_train = scoreCalc.toSklearnData(get_preds=False)
		print data_train.shape

		this_outDir = outDir + "IEX_%i_BEADS_%i_%i" % (i,j,k)

		CS.predictInteractions(scoreCalc, this_outDir, rf, num_cores, scoreF=all_scores_F,
							   verbose=True, fs=scores_to_keep)
		predF = "%s.pred.txt" % (this_outDir)
		clustering_CMD = "java -jar src/cluster_one-1.0.jar %s > %s.clust.txt" % (predF, this_outDir )
		os.system(clustering_CMD)
		eval_line = CS.clustering_evaluation([training_p, training_n, go_complexes, corum_complexes], scoreCalc, this_outDir,
								 ",".join([score.name for score in this_fs]), num_cores,
								 rf)
		out_lines.append("%i\t%i\t%s" % (i,j,eval_line))

	outFH = open(outDir + "IEX_%i_BEADS_%i.comb.eval.txt" % (i,j), "w")
	print >> outFH, "\n".join(out_lines)
	outFH.close()


def lineCount(filename):
	f = open(filename, "r")
	i = 0
	for l in f:
		i += 1
	f.close()
	return i

def benchmark():
	dir,   trainingF, holdoutF, allF, outdir = sys.argv[1:]

	lst = list(itertools.product([0, 1], repeat=8))
	clfs = ["SVM", "RF"]
	datasets = ["MSB", "SEQ", "MaxQMS1", "MaxQMS2"]

	holdout_comp = GS.Clusters(False)
	holdout_comp.read_file(holdoutF)

	train_comp = GS.Clusters(False)
	train_comp.read_file(trainingF)

	all_comp = GS.Clusters(False)
	all_comp.read_file(allF)
	scores = [CS.MutualInformation(2), CS.Bayes(3), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5),
			  CS.Pearson(), CS.Apex()]

	for clf in clfs:
		for ds in datasets:
			head = ""
			out_data = ""
			for fc in lst:
				fc = "".join(map(str,fc))
				if fc == "00000000": continue
				this_scores = []
				for i, feature_selection in enumerate(fc):
					if feature_selection == "1": this_scores.append(scores[i].name)
				this_scores = ",".join(this_scores)

				curdir = "%s%s.%s.%s" % (dir, fc, clf, ds)
				if os.path.isfile(curdir + ".train.pred.txt") == False:
					print "no file for %s.train.pred.txt" % curdir
					continue
				if lineCount(curdir + ".train.pred.txt") == 0:
					print "no predicted ppis for %s.train.pred.txt" % curdir
					continue

				if os.path.isfile(curdir + ".all.pred.txt") == False:
					print "no file for %s.all.pred.txt" % curdir
					continue
				if lineCount(curdir + ".all.pred.txt") == 0:
					print "no predicted ppis for %s.all.pred.txt" % curdir
					continue


				train_num_ppis = CS.lineCount(curdir + ".train.pred.txt")
				train_num_comp = CS.lineCount(curdir + ".train.clust.txt")
				all_num_ppis = CS.lineCount(curdir + ".all.pred.txt")
				all_num_comp = CS.lineCount(curdir + ".all.clust.txt")

				train_line, train_head = CS.clustering_evaluation(train_comp, "Train", curdir + ".train")
				holdout_line, holdout_head = CS.clustering_evaluation(holdout_comp, "Holdout", curdir + ".train")
				all_line, all_head = CS.clustering_evaluation(all_comp, "All", curdir + ".all")
				if head == "": head = "Feature set\tTrain Pred PPI\tTrain Pred Clusters\t%s\t%s\tAll Pred PPI\tAll Pred Clusters\t%s" % (train_head, holdout_head, all_head)
				out_data +=  "\n%s\t%i\t%i\t%s\t%i\t%i\t%s" % (this_scores, train_num_ppis, train_num_comp, "\t".join([train_line, holdout_line]), all_num_ppis, all_num_comp, all_line)
			outF = open("%s%s.%s.eval.txt" % ( outdir, clf, ds), "w")
			print >> outF, head + out_data
			outF.close()

"""


def Goldstandard_from_cluster_File(gsF):
		clusters = GS.Clusters(need_to_be_mapped=False)
		clusters.read_file(gsF)
		gs = GS.Goldstandard_from_Complexes("All")
		gs.complexes = clusters
		gs.make_pos_neg_ppis()
		return gs

class feature_selector:
	def __init__(self, feature_names, scoreCalc):
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

	def filter_scoreCalc(self, scoreCalc):
		filtered_scoreCalc = copy.deepcopy(scoreCalc)
		filtered_scoreCalc.header = np.array(filtered_scoreCalc.header)[self.to_keep_header]
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[:, self.to_keep_score]
		prot_indices = [None]* len(scoreCalc.ppiToIndex.keys())
		for p, i in scoreCalc.ppiToIndex.items(): prot_indices[i] = p
		val_rows = np.where(np.apply_along_axis( self.valid_score, 1, filtered_scoreCalc.scores)==True)[0]
		filtered_scoreCalc.scores = filtered_scoreCalc.scores[val_rows, :]
		prot_indices = np.array(prot_indices)[val_rows]
		filtered_scoreCalc.ppiToIndex = {}
		filtered_scoreCalc.IndexToPpi = {}
		for i, p in enumerate(prot_indices):
			filtered_scoreCalc.ppiToIndex[p] = i
			filtered_scoreCalc.IndexToPpi[i] = p
 		return filtered_scoreCalc

	def get_next(self):
		edge, scores = self.scoreCalc.get_next()
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
#					,[CS.Bayes(1), CS.Bayes(2)]
#					,[CS.Bayes(1), CS.Bayes(3)]
#					,[CS.Bayes(2), CS.Bayes(3)]
#					,[CS.Bayes(1), CS.Bayes(2), CS.Bayes(3)]]

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
		all_gs = Goldstandard_from_cluster_File(ref_complexes)

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
	feature_comb = feature_selector([fs.name for fs in feature_combination], scoreCalc)
	print feature_comb.scoreCalc.scores.shape

	print "Num valid filterd ppis: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()))
	print "Num valid filtered pos: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()) & set(all_gs.positive))
	print "Num valid filtered negative: %i" % len(set(feature_comb.scoreCalc.ppiToIndex.keys()) & set(all_gs.negative))

	out_prefix = "_".join([fs.name for fs in feature_combination])
	train, eval = all_gs.split_into_holdout_training(set(feature_comb.scoreCalc.ppiToIndex.keys()), no_overlapp=no_overlap_in_training)

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
	scores = [CS.MutualInformation(2), CS.Bayes(2), CS.Euclidiean(), CS.Wcc(), CS.Jaccard(), CS.Poisson(5), CS.Pearson(), CS.Apex()]
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


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass

	#11000100 (MI, Bayes, PCC+N)