#!/usr/bin/env python
from __future__ import division
import sys
import os
import re
import subprocess
import gzip
import math
import numpy as np
import CalculateCoElutionScores as CS
import copy

def getSubSet(eData, indices):
	eData.elutionMat = np.array(eData.elutionMat)
	eData.elutionMat = eData.elutionMat[:, indices]
	eData.normedElutionMat = eData.normedElutionMat[:, indices]
	not_empty_rows = np.where(eData.elutionMat.any(1))[0]
#	eData.elutionMat = eData.elutionMat[not_empty_rows, :] do not remove empty rows else index two prot needs to be readjusted
	prots_todel = []
	for prot in eData.prot2Index:
		prot_index = eData.prot2Index[prot]
		if prot_index not in not_empty_rows: prots_todel.append(prot)

	for prot in prots_todel:
		del eData.prot2Index[prot]

	return eData

def get_sim_eData(eD1, eD2, score, cut_off=0.5):
	sharedprots = set(eD1.prot2Index.keys()) & set(eD2.prot2Index.keys())
	matched_prots = set([])
	scores = []
	score1 = copy.deepcopy(score)
	score2 = copy.deepcopy(score)
	for prot in sharedprots:
		profileD1, _ = score1.getScores(prot, prot, eD1)
		score.clear()
		profileD2, _ = score2.getScores(prot, prot, eD2)
		score.clear()
		sim_score = score.calculateScore(profileD1, profileD2)
		scores.append(sim_score)
	return sum(scores)/len(scores)
#		if sim_score > cut_off:
#			matched_prots.add(prot)
#	return len(matched_prots)/len(sharedprots)

def main():
	experimentF, outDir = sys.argv[1:]

	experiments = CS.ElutionData(experimentF)

	ref_exp = getSubSet(copy.deepcopy(experiments), range(10))
	tmt1 = getSubSet(copy.deepcopy(experiments), range(10, 20))
	tmt2 = getSubSet(copy.deepcopy(experiments), range(20, 30))
	tmt3 = getSubSet(copy.deepcopy(experiments), range(30, 40))
	tmt_combined = getSubSet(copy.deepcopy(experiments), range(10,40))

	get_sim_eData(ref_exp, tmt1, CS.Apex())

	experiments = [ref_exp, tmt1, tmt2, tmt3]

	exp_names = ["Ref", "TMT1", "TMT2", "TMT3"]

	simMat = [[0]*len(experiments)]*len(experiments)

	scores = [CS.Apex(), CS.Wcc(), CS.MutualInformation(), CS.Euclidiean(), CS.Pearson(), CS.Jaccard()]
#	scores = [CS.Wcc(), CS.Euclidiean()]


	for st in scores:
		datFH = open(outDir + ".%s" % st.name + ".txt", "w")
		datFH.write("\t"+"\t".join(exp_names) + "\n")
		for i in range(len(experiments)):
			datFH.write(exp_names[i])
			for j in range(len(experiments)):
				if i == j:
					datFH.write("\t1.0")
				else:
					datFH.write("\t" + str(get_sim_eData(experiments[i], experiments[j], st)))
			datFH.write("\n")
		datFH.close()
		cmd =  "src/plotHeatmap.R %s.%s.txt %s %s.%s.pdf" % (outDir, st.name, st.name, outDir, st.name)
		os.system(cmd)


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
