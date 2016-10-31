#!/usr/bin/env python
from __future__ import division
import CalculateCoElutionScores as CS
import GoldStandard as GS
import copy
import numpy as np
import sys
import os
from matplotlib import pyplot as plt
#from matplotlib_venn import venn3


def mapPPI(ppiF, mappingF, outF, foundprots=""):
	def readMapping(mappingF):
		mappingFH = open(mappingF)
		out = {}
		todel = set([])
		mappingFH.readline()
		for line in mappingFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			if len(linesplit) <= 1 : continue
			ida, idbs = line.rstrip().split("\t")
			if foundprots != "" and ida not in foundprots: continue
			idbs = idbs.split(";")
			for idb in idbs:
				if idb=="":continue
				if idb not in out:
					out[idb] = ida
				else:
					todel.add(idb)


		for id in todel:
			del out[id]

		return out

	if mappingF !="NA":
		mapping = readMapping(mappingF)
	else:
		mapping = foundprots

	outFH = open(outF, "w")
	ppiFH = open(ppiF)
	for line in ppiFH:
		line = line.rstrip()
		linesplit = line.split("\t")
		ida, idb = linesplit[0:2]
		if ida == idb: continue
		if mapping!="" and (ida not in mapping or idb not in mapping): continue
		edge = "\t".join(sorted([ida, idb]))
		if  mappingF !="NA":
			ida_mapped = mapping[ida]
			idb_mapped = mapping[idb]
			edge =  "\t".join(sorted([ida_mapped, idb_mapped]))
		print >> outFH, "%s\t%s" % (edge, "\t".join(linesplit[2:]))
	outFH.close()
	ppiFH.close()





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

def test_geneMania():
	scorecalc = CS.CalculateCoElutionScores()
	scorecalc.readTable("/Users/florian/workspace/scratch/EPIC_out/tmp2")
	ppis =  scorecalc.ppiToIndex.keys()
	target_species = "6239"
	genemania = CS.Genemania(target_species)

def readPPIs(ppiF):
	ppiFH = open(ppiF)
	out = set()
	for line in ppiFH:
		out.add(line.rstrip())
	return out

def main():
	test_geneMania()
	sys.exit()

	elutionFH = open("/Users/florian/workspace/scratch/EPIC_out/1D_files.txt")
	elutionProts = set([])
	for elutionFile in elutionFH:
		elutionFile = elutionFile.rstrip()
		elutionData = CS.ElutionData(elutionFile)
		elutionProts = elutionProts | set(elutionData.prot2Index.keys())



	#get_cluster_overlapp()
	gs_cluster_overlapp()
	sys.exit()

#	ppiF, mappingF, outF = sys.argv[1:]
#	mapPPI(ppiF, mappingF, outF, elutionProts)
#	sys.exit()

	netA, netB, netC = sys.argv[1:]
	netA = readPPIs(netA)
	netB = readPPIs(netB)
	netC = readPPIs(netC)

	net_overlap = [len(netA), len(netB), len(netA&netB), len(netC), len(netA&netC), len(netB&netC), len(netA&netB&netC)]
	venn3([netA, netB, netC], ("STRING", "IntAct", "BIOGRID"))
	plt.savefig("/Users/florian/workspace/scratch/EPIC_out/GS_net.pdf", bbox_inches="tight")
	plt.close()




if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass
