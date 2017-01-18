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
#from matplotlib_venn import venn3
import glob
import random as re
import fs as fs



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

def epic_read_eval(fs_eval_F):
	def norm_score(scores, i):
		min_score = min(scores[:,i])
		max_score = max(scores[:,i])
		for j in range(len(scores[:,i])):
			scores[j, i]  = (scores[j,i]-min_score)/(max_score-min_score)
	fs_eval_FH = open(fs_eval_F)
	vals = []
	fs = []
	header = fs_eval_FH.readline().strip().split("\t")[7:]

	for line in fs_eval_FH:
		line = line.rstrip().split("\t")
		vals.append(map(float, line[7:]))
		fs.append(line[0])
	fs_eval_FH.close()

	vals = np.array(vals)
	zvals = zscore(vals)
#	zvals = copy.copy(vals)
#	norm_score(zvals, 0)
#	norm_score(zvals, 1)

	return header, fs, vals, zvals

def EPICFS_eval(fs_eval_F="", verbose=True, filter=True):
	if fs_eval_F =="": fs_eval_F = sys.argv[1]

	header, fs, vals = epic_read_eval(fs_eval_F)

	if filter:
		sel_vals = list(np.where(vals[:, 1] >= 100)[0])
		vals = vals[sel_vals,]

	best_fs = []

	for i in range(len(header)):
		max_index = np.argmax(np.array(vals[:,i]))
		pred_comp_val = str(int(vals[max_index, 1]))
		best_fs.append(fs[max_index])

		if not verbose:continue
		if i <=1:
			print "%i (%s)\t%s" % (int(vals[max_index, i]), pred_comp_val, fs[max_index])
		else:
			print "%.3f (%s)\t%s" % (vals[max_index,i], pred_comp_val, fs[max_index])

	return best_fs

def get_max_vals(vals, fs, header):
	out = {}
	for i in range(len(header)):
		max_index = np.argmax(np.array(vals[:,i]))
		out[header[i]] =  vals[max_index,i]
	return out

def EPIC_count_fs():

	counts = {}
	best_fs = EPICFS_eval(verbose=False)

	for cat in [[2],[3,4,5],[6],[7],[8],[9],[10],[11,12,13],[14],[15],[16],[17]]:
		fs = set()
		for i in cat:
			fs |= set(best_fs[i].split(","))
		for feature in fs:
			if feature not in counts:counts[feature]=0
			counts[feature]+= 1

	for feature in sorted(counts.keys()):
		print feature

	for feature in sorted(counts.keys()):
		print counts[feature]

def EPIC_eval_fs_DIST():
	def getScore(scores):
		out = []
		for cat in [[0, 1], [2,9], [3, 4, 5], [6, 7, 8]]:
			out.append(sum(scores[cat])/len(cat))
		return out
	def dist(a,b):
		return distance.euclidean(a,b)

	all_scores = {}

	fs_eval_Files = sys.argv[1:]
	max_vals = []
	for fs_eval_F in fs_eval_Files:
		header, fs, vals, zvals = epic_read_eval(fs_eval_F)
		fs = np.array(fs)
#		sel_vals = list(np.where(vals[:, 1] >= 100)[0])
#		vals = vals[sel_vals,]
#		zvals = zvals[sel_vals,]
#		fs = np.array(fs)[sel_vals,]
		all_scores[fs_eval_F] = (header, fs, vals, zvals)



		this_max_vals = []
		for i in range(len(header)):
			max_index = np.argmax(np.array(zvals[:, i]))
			this_max_vals.append(zvals[max_index, i])
		max_vals.append(getScore(np.array(this_max_vals)))

	max_vals = np.array(max_vals)

	tmp = []
	for i in range(max_vals.shape[1]):
		tmp.append(max(max_vals[:,i]))
	max_vals = np.array(tmp)

	scores = {}
	for score in all_scores:
		header, fs, vals, zvals = all_scores[score]
		for i in range(len(fs)):
			this_f = fs[i]
			this_vals = getScore(zvals[i,:])
			this_dist = dist(this_vals,max_vals)
			scores[(score, this_f)] = this_dist

	best_dist = min(scores.values())
	for score, f in scores:
		header, fs, vals, zvals = all_scores[score]
		if scores[(score,f)]==best_dist:
			print score
			print f
			print "\n".join(map(str, vals[np.where(fs == f)[0], :][0]))

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
			vals.append(map(float, line[7:]))
			fs.append(line[0])
		fs_eval_FH.close()
	vals = np.array(vals)
	for row in np.corrcoef(np.transpose(vals)):
		print "\t".join(map(str, row))


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

def bac_EPIC():
	elutionFile, inparanoid_file, cluster_file, cluster_map_file = sys.argv[1:]

	eData = CS.ElutionData(elutionFile)
	foundprots = set(eData.prot2Index.keys())
	orthmapper = GS.Inparanoid(taxid="NA", inparanoid_cutoff=0.95, foundProts=foundprots)
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

	gs = GS.Goldstandard_from_Complexes("Intact")
	gs.make_reference_data(GS.Intact_clusters(species="escherichia_coli"), orthmapper, found_prots=foundprots)

	gs = GS.Goldstandard_from_Complexes("GO")
	gs.make_reference_data(GS.QuickGO(taxid="83333"), orthmapper, found_prots=foundprots)



def main():
#	EPIC_cor()
#	EPICFS_eval()
#	EPIC_count_fs()
#	EPIC_eval_fs()
#	EPIC_eval_fs_DIST()
#	exp_comb()
#	bac_EPIC()

	foundprots, elution_datas = CS.load_data("/Users/florian//workspace/scratch/EPIC_out/input/elution_profiles/MSB", [])

	_, _, _, _, go_complexes, corum_complexes, intact_complexes	= CS.create_goldstandard("6239",foundprots)
	outFH = open("/Users/florian/workspace/scratch/EPIC_out/go.comp.txt", "w")
	print >> outFH, go_complexes.to_string()
	outFH.close()

	outFH = open("/Users/florian/workspace/scratch/EPIC_out/corum.comp.txt", "w")
	print >> outFH, corum_complexes.to_string()
	outFH.close()


	outFH = open("/Users/florian/workspace/scratch/EPIC_out/intact.comp.txt", "w")
	print >> outFH, intact_complexes.to_string()
	outFH.close()


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
