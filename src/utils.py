#!/usr/bin/env python
from __future__ import division

from scipy.spatial import distance
from scipy.stats import zscore

import CalculateCoElutionScores as CS
import GoldStandard as GS
import copy
import numpy as np
import sys
import os, math
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import glob
import random as re
import fs as fs
import itertools


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

def reza_overlap():
	dataF = "/Users/florian/workspace/mouseMap/data/MOUSE/proteinGroups_BRAIN_SC_modified.txt"
	tmt = ["HS", "LS", "161115_L", "161115_S"]
	for j,i in enumerate([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]):
		dataFH = open(dataF)
		dataFH.readline()
		prots = [set([]), set([]), set([])]
		for line in dataFH:
			line = np.array(line.rstrip().split("\t"))
			id = line[0]
			cyt, mem, neu = line[i]
			if float(cyt) > 1: prots[0].add(id)
			if float(mem) > 1: prots[1].add(id)
			if float(neu) > 1: prots[2].add(id)
		dataFH.close()
		venn3(prots, ("Cytosolic", "Membrane", "Nucleus"))
		plt.title(tmt[j], fontsize=14)
		plt.savefig("/Users/florian/workspace/scratch/%s.pdf" % tmt[j], bbox_inches="tight")
		plt.close()







def gs_cluster_overlapp():
	elutionFH = open("/Users/florian/workspace/scratch/EPIC_out/1D_files.txt")
	elutionProts = set([])
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = CS.ElutionData(elutionFile)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())


	target_taxid = "9606"

	corum_complexes = GS.CORUM()
	corum_gs = GS.Goldstandard_from_Complexes("CORUM")
	corum_gs.make_reference_data(corum_complexes, target_taxid, "")
	outFH = open("/Users/florian/workspace/scratch/EPIC_out/CORUM.txt", "w")
	print >> outFH, corum_gs.complexes.to_string()

	intact_complexes = GS.Intact_clusters()
	intact_gs = GS.Goldstandard_from_Complexes("Intact")
	intact_gs.make_reference_data(intact_complexes, target_taxid, "")

	outFH = open("/Users/florian/workspace/scratch/EPIC_out/Intact.txt", "w")
	print >> outFH, intact_gs.complexes.to_string()

	go_complexes = GS.QuickGO(target_taxid)
	go_gs = GS.Goldstandard_from_Complexes("GO")
	go_gs.make_reference_data(go_complexes, target_taxid, "")
	outFH = open("/Users/florian/workspace/scratch/EPIC_out/GO.txt", "w")
	print >> outFH, go_gs.complexes.to_string()
	sys.exit()

	corum_p, corum_n = corum_gs.get_goldstandard()
	intact_p, intact_n = intact_gs.get_goldstandard()
	go_p, go_n = go_gs.get_goldstandard()

	prot_overlap = [len(go_prots),len(corum_prots),len(go_prots&corum_prots),len(intact_prots),len(go_prots & intact_prots),len(corum_prots&intact_prots),len(go_prots&intact_prots&corum_prots)]
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_prots_raw.pdf")
	plt.close()


	go_prots = set(go_gs.get_complexes().getProtToComplexMap().keys())
	intact_prots = set(intact_gs.get_complexes().getProtToComplexMap().keys())
	corum_prots =  set(corum_gs.get_complexes().getProtToComplexMap().keys())

	venn3([go_prots, corum_prots, intact_prots], ("GO", "CORUM", "Intact"))
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_prots.pdf")
	plt.close()

	venn3([go_p, corum_p, intact_p] , ("GO", "CORUM", "Intact"))
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_p.pdf", bbox_inches="tight")
	plt.close()

	venn3([go_n, corum_n, intact_n] , ("GO", "CORUM", "Intact"))
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_n.pdf", bbox_inches="tight")
	plt.close()

	outline= "\tGO\tCORUM\tIntAct\n"
	for compA in [go_gs.get_complexes(), corum_gs.get_complexes(), intact_gs.get_complexes()]:
		for compB in [go_gs.get_complexes(), corum_gs.get_complexes(), intact_gs.get_complexes()]:
			outline += "\t" + str(compA.getOverlapp(compB))
		outline += "\n"

	print outline

	scorecalc = CS.CalculateCoElutionScores()
	scorecalc.readTable("/Users/florian/workspace/scratch/EPIC_out/All.scores.sel.txt")
	found_ppis = set(scorecalc.ppiToIndex.keys())
	go_p = found_ppis & go_p
	corum_p = found_ppis & corum_p
	intact_p = found_ppis & intact_p

	venn3([go_p, corum_p, intact_p] , ("GO", "CORUM", "Intact"))
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_p_val.pdf")
	plt.close()

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
#	zvals = copy.copy(vals)
#	norm_score(zvals, [0,1])
#	norm_score(zvals, range(10))

	return header, fs, vals, zvals

def EPIC_eval_fs_DIST():

	def getScore(scores):
		out = []
#		for cat in [ [1], [2],[3],[4],[5], [6], [7], [9], [9]]:
#		for cat in [[1], [2, 9], [3, 4], [6], [7]]:
#		for cat in [[1], [2], [9], [3, 4, 5], [6, 7, 8]]:
#		for cat in [[1,2,9], [3,4],[6],[7]]:
#		for cat in [[0], [1, 2, 9],[ 3, 4, 5, [6], [7]]:
		for cat in [[1, 2, 9], [3,4,5,6], [7]]:
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

def EPIC_eval_fs():
	fs_a, fs_b, fs_eval_F = sys.argv[1:]

	fs_a = fs_a.split("_")
	fs_a = [set(f.split(",")) for f in fs_a]

	fs_b = fs_b.split("_")
	fs_b = [set(f.split(",")) for f in fs_b]
	header, fs, vals = epic_read_eval(fs_eval_F)
	sel_vals = list(np.where(vals[:, 1] >= 100)[0])
	vals = vals[sel_vals,]
	fs = np.array(fs)[sel_vals,]

	fs_index = {}
	for i in range(len(fs)):
		fs_index[tuple(set(fs[i].split(",")))] = i

	best = {}
	for h in header[2:]:
		best[h] = (set([]), 0.0)
	cur_fs = set([])
	for f in fs_a:
		cur_fs |= f
		tmp_fs = tuple(set([]) | cur_fs)
		if tmp_fs not in fs_index: continue
		print cur_fs
		cur_vals = vals[fs_index[tmp_fs],:]
		for i, val in enumerate(cur_vals):
			if i <2:continue
			cat = header[i]
			if val > best[cat][1]:
				best[cat] = (tmp_fs, val)


	cur_fs = set([])
	for f in fs_b:
		cur_fs |= set(f)
		tmp_fs = tuple(set([]) | cur_fs)
		if tmp_fs not in fs_index: continue
		print cur_fs
		cur_vals = vals[fs_index[tmp_fs], :]
		for i, val in enumerate(cur_vals):
			if i < 2: continue
			cat = header[i]
			if val > best[cat][1]:
				best[cat] = (tmp_fs, val)

	counts = {}
	for cat in [[2],[3,4,5],[6],[7],[8],[9],[10],[11,12,13],[14],[15],[16],[17]]:
		cats = [header[i] for i in cat]
		tmp_counts={}
		for cat in cats:
			f, val = best[cat]
			if tuple(f) not in tmp_counts:tmp_counts[tuple(f)]=0
			tmp_counts[tuple(f)]+=1

		max_val = max(tmp_counts.values())
		for f in tmp_counts:
			this_v = tmp_counts[f]
			if max_val == this_v:
				print f
				print cats
				if f not in counts:counts[f]=0
				counts[f]+=1

	max_val = max(counts.values())
	for f in counts:
		if counts[f] == max_val:
			print ",".join(f)
			print counts[f]
			print "\n".join(map(str,vals[fs_index[f], :]))


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

def count_enrichments():
	dir = "/Users/florian/workspace/HuRI/output/"
	files = ["cofraccomplexes.go_enrichment.txt", "bioplexcomplexes.go_enrichment.txt", "C1_CoFrac_clust.go_enrichment.txt", "C1_Bioplex_clust.go_enrichment.txt", "C1_UNION_QBC_clust.go_enrichment.txt"]
	names = ["cofrac", "bioplex", "C1_cofrac", "C1_bioplex", "C1_union"]
	labels = {"BP": 1, "MF": 2, "CC": 3}
	outF = sys.argv[1]
	out = []
	for i, file in enumerate(files):
		out_c = np.array([0, 0, 0,0])
		all_counts = []
		name = names[i]
		num_comp = CS.lineCount(dir+file)
		inFH = open(dir+file)
		for line in inFH:
			line = line.split("\t")
			size = len(line[0].split(";"))
			if size < 3: continue
#			if size > 7: size = "7+"
#			size = str(size)
			go_counts = []

			for label in labels:
				index = labels[label]
				counts = len(line[index].split(";"))
				out.append("%s\t%s\t%f" % (name, label, counts/size))
				out_c[index] += counts
		print(file + "\t" + "\t".join(map(str,out_c[1:]/num_comp)))
		inFH.close()

	outFH = open(outF, 'w')
	print >> outFH, "Clusters\tGO category\tvalue"
	print >> outFH, "\n".join(out)
	outFH.close()

def mk_overlapGraph():
	clustF, outDir = sys.argv[1:]

	def simco(a, b):
		tmpa = set(a)
		tmpb = set(b)
		return len(tmpa & tmpb) / (min(len(tmpa), len(tmpb)))


	def overlap(a, b):
		tmpa = set(a)
		tmpb = set(b)
		overlap = math.pow((len(tmpa & tmpb)), 2) / (len(tmpa) * len(tmpb))
		return overlap

	comp = GS.Clusters(False)
	comp.read_file(clustF)

	simco_outF = open(outDir + ".simco.edges.txt", "w")
	overlap_outF = open(outDir + ".overlap.edges.txt", "w")

	for i, comp_i in comp.complexes.items():
		for j, comp_j in  comp.complexes.items():
			if j <= i: continue
			this_simco =  simco(comp_i, comp_j)
			if this_simco > 0:
				print >> simco_outF, "%i\t%i\t%f" % (i, j, this_simco)
			this_overlap = overlap(comp_i, comp_j)
			if this_simco > 0:
				print >> overlap_outF, "%i\t%i\t%f" % (i, j, this_overlap)

	simco_outF.close()
	overlap_outF.close()


def orthmap_complexes():
	compF, mapF, outF = sys.argv[1:]
	comps = GS.Clusters(False)
	comps.read_file(compF)
	print len(comps.complexes)

	cluster_map = {}
	cluster_mapFH = open(mapF)
	for line in cluster_mapFH:
		line = line.rstrip()
		protid, compid = line.split("\t")
		cluster_map[compid] = protid
	cluster_mapFH.close()

	for cluster in comps.complexes:
		new_ids = set([])
		for cluster_id in  comps.complexes[cluster]:
			if cluster_id in cluster_map:
				protid = cluster_map[cluster_id]
				new_ids.add(protid)
			else:
				print "No map for %s" % cluster_id
		comps.complexes[cluster] = new_ids
	outFH = open(outF, "w")
	print >> outFH, comps.to_string()
	outFH.close()



def bac_EPIC():
	elutionFile, inparanoid_file, cluster_file, cluster_map_file, outDir = sys.argv[1:]

	scores = [CS.Wcc(), CS.MutualInformation(), CS.Poisson(5), CS.Jaccard(), CS.Apex()]
	eData = CS.ElutionData(elutionFile)
	for s in scores: s.init(eData)
	foundprots = set(eData.prot2Index.keys())
	orthmapper = GS.Inparanoid(taxid="NA", inparanoid_cutoff=0.95, foundProts="")
	orthmapper.readTable(inparanoid_file)
	ecocyc_compl = GS.Clusters(need_to_be_mapped=True)
	ecocyc_compl.read_file(cluster_file)
	ecocyc_compl.lb=2

	cluster_map = {}
	cluster_mapFH = open(cluster_map_file)
	for line in cluster_mapFH:
			line = line.rstrip()
			compid, protid = line.split("\t")
			cluster_map[compid] = protid
	cluster_mapFH.close()

	for cluster in ecocyc_compl.complexes:
		for cluster_id in  ecocyc_compl.complexes[cluster]:
			if cluster_id in cluster_map:
				protid = cluster_map[cluster_id]
				ecocyc_compl.complexes[cluster].remove(cluster_id)
				ecocyc_compl.complexes[cluster].add(protid)
	print "Raw: %i" % len(ecocyc_compl.complexes.keys())

	ecocyc_compl.filter_complexes()
	print "Filter: %i" % len(ecocyc_compl.complexes.keys())

	ecocyc_compl.merge_complexes()
	ecocyc_compl.filter_complexes()
	print "Merge: %i" % len(ecocyc_compl.complexes.keys())

	orthmapper.mapComplexes(ecocyc_compl)
	ecocyc_compl.filter_complexes()
	print "Orthmap: %i" % len(ecocyc_compl.complexes.keys())

	ecocyc_compl.remove_proteins(foundprots)
	ecocyc_compl.filter_complexes()
	print "Foundprots: %i" % len(ecocyc_compl.complexes.keys())

	gs_p, gs_n = ecocyc_compl.getPositiveAndNegativeInteractions()

	print len(gs_p)
	print len(gs_n)
	print gs_p

	scoreCalc = CS.CalculateCoElutionScores()
	scoreCalc.positive = gs_p
	scoreCalc.negative = gs_n
	scoreCalc.readTable(outDir + ".scores.txt")
	scoreCalc.rebalance()
	print scoreCalc.scores.shape
#	scoreCalc.calculate_coelutionDatas([eData], scores, outDir, 4)
	CS.bench_scores(scoreCalc, outDir, 4, useForest=True)

	CS.predictInteractions(scoreCalc, outDir, True, 4, scoreF="", verbose=True, fs="", score_cutoff=0)
	predF = "%s.pred.txt" % (outDir)
	clustering_CMD = "java -jar cluster_one-1.0.jar %s > %s.clust.txt" % (predF, outDir)
	print(clustering_CMD)
	os.system(clustering_CMD)

#	gs = GS.Goldstandard_from_Complexes("Intact")
#	gs.make_reference_data([GS.Intact_clusters(species="escherichia_coli"), GS.QuickGO(taxid="83333")], orthmapper, found_prots=foundprots)

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


def main():
#	EPIC_cor()
#	EPIC_count_fs()
#	EPIC_eval_fs()
#	EPIC_eval_fs_DIST()
#	exp_comb()
#	bac_EPIC()
#	orthmap_complexes()
#	count_enrichments()
# 	reza_overlap()
#	benchmark()
#	mk_overlapGraph()
	make_ref("/Users/florian/workspace/scratch/EPIC_out/input/elution_profiles/MSB/","6239")

def make_ref(e_file_dir, taxid):
	foundprots, elution_datas = CS.load_data(e_file_dir, [])
	CS.create_goldstandard(taxid, foundprots)



"""
	for name, comp_ref in [("corum", GS.CORUM()), ("intact", GS.Intact_clusters()), ("go", GS.QuickGO("9606"))]:
		comp_ref.need_to_be_mapped = False
		outFH = open("/Users/florian/workspace/scratch/EPIC_out/%s.raw.comp.txt" % name, "w")
		print >> outFH, comp_ref.get_complexes().to_string()
		outFH.close()
		gs =  GS.Goldstandard_from_Complexes(name)
		gs.make_reference_data([comp_ref], "", found_prots="")
		outFH = open("/Users/florian/workspace/scratch/EPIC_out/%s.processed.comp.txt" % name, "w")
		print >> outFH, gs.complexes.to_string()
		outFH.close()

	_, _, all_p, all_n, _, _, all_comp	= CS.create_goldstandard("9606","")


	outFH = open("/Users/florian/workspace/scratch/EPIC_out/combined.all_p.txt", "w")
	print >> outFH, "\n".join(all_p)
	outFH.close()

	outFH = open("/Users/florian/workspace/scratch/EPIC_out/combined.all_n.txt", "w")
	print >> outFH, "\n".join(all_n)
	outFH.close()

	outFH = open("/Users/florian/workspace/scratch/EPIC_out/combined.comp.txt", "w")
	print >> outFH, all_comp.to_string()
	outFH.close()
"""


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
