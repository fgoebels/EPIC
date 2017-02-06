#!/usr/bin/python
from __future__ import division
import numpy as np
import urllib2
from xml.dom import minidom
import os
import copy
import re
import math
import sys


class Goldstandard_from_PPIs():
	def __init__(self, refF):
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
		for i in range(len(allprots)):
			prot_i = allprots[i]
			for j in range(i+1, len(allprots)):
				prot_j = allprots[j]
				edge = "\t".join(sorted([prot_i, prot_j]))
				if edge in self.goldstandard_positive: continue
				self.goldstandard_negative.add(edge)


class Goldstandard_from_reference_File():
	def __init__(self, refF, found_prots=""):
		self.goldstandard_positive, self.goldstandard_negative = self.read_reference_file(refF, found_prots)

	def read_reference_file(self, refF, found_prots=""):
			pos = set([])
			neg = set([])
			refFH = open(refF)
			for line in refFH:
					line = line.rstrip()
					(idA, idB, label) = line.split("\t")
					if found_prots != "":
							if idA not in found_prots: continue
							if idB not in found_prots: continue
					edge = "\t".join(sorted([idA, idB]))
					if label == "negative": neg.add(edge)
					if label == "positive": pos.add(edge)
			return pos, neg

class Goldstandard_from_cluster_File():
	def __init__(self, gsF, found_prots=""):
		self.clusters = ""
		self.goldstandard_positive, self.goldstandard_negative = self.readGS(gsF, found_prots)

	def readGS(self, gsF, found_prots):
		clusters = Clusters(need_to_be_mapped=False)
		clusters.read_file(gsF)
		clusters.filter_complexes()
		clusters.merge_complexes()
		clusters.filter_complexes()
		if found_prots !="": clusters.remove_proteins(found_prots)
		self.clusters = clusters
		clusters.filter_complexes()
		self.clusters = clusters
		return clusters.getPositiveAndNegativeInteractions()


class Goldstandard_from_Complexes():

	def __init__(self, name="unnamed"):
		self.complexes = ""
		self.name = name
		self.goldstandard_positive, self.goldstandard_negative = set([]), set([])

	def make_reference_data(self, db_clusters, orthmap="", found_prots=""):

		self.complexes = Clusters(False)
		total_complexes = 0

		for db_clust in db_clusters:
			tmp_clust = copy.deepcopy(db_clust.get_complexes())
			total_complexes += len(tmp_clust.complexes.keys())
			if tmp_clust.need_to_be_mapped == True and orthmap !="":
				orthmap.mapComplexes(tmp_clust)
			for compl in tmp_clust.complexes:
				self.complexes.addComplex(compl, tmp_clust.complexes[compl])

		print "Total number of complexes %i in %s" % (total_complexes, self.name)
		print "Number of complexes after ortholog mapping %i complexes in %s" % (len(self.complexes.complexes), self.name)


		self.complexes.filter_complexes()
		print "After size filtering %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		self.complexes.merge_complexes()
		self.complexes.filter_complexes()
		print "After mergning %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		if found_prots != "":
			self.complexes.remove_proteins(found_prots)
			self.complexes.lb = 2
			self.complexes.filter_complexes()
			print "After removing not indetified proteins %i number of complexes in % s" % (len(self.complexes.complexes), self.name)

		self.goldstandard_positive, self.goldstandard_negative = self.complexes.getPositiveAndNegativeInteractions()

	def get_complexes(self):
		return self.complexes

	def get_goldstandard(self):
		return self.goldstandard_positive, self.goldstandard_negative

class Intact_clusters():

	def __init__(self, species = "homo_sapiens"):
		self.complexes = Clusters(need_to_be_mapped=True)
		self.need_to_be_mapped = True
		self.load_data(species)


	def get_complexes(self):
		return self.complexes

	def load_data(self, species):
		intact_url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/%s.tsv" % species
		intact_url_FH = urllib2.urlopen(intact_url)
		intact_url_FH.readline()
		i = 0
		for line in intact_url_FH:
			line = line.rstrip()
			linesplit = line.split("\t")
			evidence = linesplit[5]
			if not evidence.startswith("ECO:0000353"): continue
			members = linesplit[4]
			members = re.sub("\(\d+\)", "", members)
			members = re.sub("-\d+", "", members)
			members = re.sub("-PRO_\d+", "", members)
			members = set(members.split("|"))
			self.complexes.addComplex(i, members)
			i += 1

# @author Florian Goebels
# Wrapper class for downloading and handling CORUM protein complex information taken from here: http://mips.helmholtz-muenchen.de/genre/proj/corum/
class CORUM():
	# @author Florian Goebels
	# object constructor
	# @param
	#		lb lower bound complex should have at least lb members
	#		ub upper bound complex should have at most  ub members
	#		overlap_cutoff merge complexes that have an overlap_score > overlap_cutoff
	#		source_species select for which species the complexes should be maintained
	def __init__(self, source_species_regex = "(Human|Mammalia)"):
		# static regex for identifying valid bochemical evidences codes
		self.biochemical_evidences_regex ="MI:(2193|2192|2191|2197|2195|2194|2199|2198|0807|0401|0400|0406|0405|0404|0089|0084|0081|0007|0006|0004|0513|1029|0979|0009|0008|0841|1312|2188|2189|0411|0412|0413|0928|0415|0417|0098|0729|0920|0921|0603|0602|0605|0604|0402|0095|0096|0606|0091|0092|1142|1145|1147|0019|1309|0696|0697|0695|0858|0698|0699|0425|0424|0420|0423|0991|0990|0993|0992|0995|0994|0997|0996|0999|0998|1028|1011|1010|1314|0027|1313|0029|0028|0227|0226|0225|0900|0901|0430|0434|0435|1008|1009|0989|1004|1005|0984|1007|1000|0983|1002|1229|1087|1325|0034|0030|0031|0972|0879|0870|1036|0678|1031|1035|1034|0676|0440|1138|1236|0049|0048|1232|0047|1137|0419|0963|1026|1003|1022|0808|0515|0514|1187|0516|0511|1183|0512|0887|0880|0889|0115|1006|1249|0982|0953|1001|0508|0509|0657|0814|1190|1191|0813|0066|0892|0899|1211|0108|1218|1352|1354|0949|0946|0947|0073|0071|1019|2168|0700|2167|1252|1017|0276|1189|1184)"
		self.source_species_regex = source_species_regex
		self.getCORUM()
		self.readCORUM()


	def get_complexes(self):
		return self.complexes

	# @author Florian Goebels
	# downloads current version of corum and safe it to wd/data folder as corum.txt
	def getCORUM(self):
		self.corum_raw = {}
		corum_url = "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt"
		corum_url_FH = urllib2.urlopen(corum_url)
		for line in corum_url_FH:
			line = line.rstrip()
			linesplit = np.array(line.split("\t"))
			corum_id = linesplit[0]
			self.corum_raw[corum_id] = linesplit[(2,5,7),]

	# @author Florian Goebels
	# reads in CORUM from flat file
	def readCORUM(self):
		self.complexes = Clusters(need_to_be_mapped=True)
		for comp in self.corum_raw:
			(species, prots, evidence) = self.corum_raw[comp]
			# bool(...) returns true if evidence code is found => not bool(...) is true if not valid evidence is found, and thus skip this complex
			if not bool(re.search(self.biochemical_evidences_regex, evidence)):
				continue
			if not bool(re.search(self.source_species_regex, species)): continue
			prots = set(prots.split(";"))
			self.complexes.addComplex(comp, prots)

class Clusters():

	def __init__(self, need_to_be_mapped, overlap_cutoff = 0.8, lb = 3, ub = 50):
		self.complexes = {}
		self.overlap_cutoff = overlap_cutoff
		self.ub = ub
		self.lb = lb
		self.need_to_be_mapped = need_to_be_mapped

	def addComplex(self, complex, members):
		if complex not in self.complexes: self.complexes[complex] = set([])
		self.complexes[complex] = self.complexes[complex] | members

	def read_file(self, clusterF):
		clusterFH = open(clusterF)
		i = 0
		for line in clusterFH:
			line = line.rstrip()
			prots = set(line.split("\t"))
			self.addComplex(i, prots)
			i+=1

	def write_cuslter_file(self, outF):
		outFH = open(outF, "w")
		print >> outFH, self.to_string()
		outFH.close()

	def to_string(self):
		out = []
		for clust in self.complexes:
			prots = self.complexes[clust]
			out.append("\t".join(prots))
		return "\n".join(out)

	# @author Florian Goebels
	# creats all possible positive and negative protein interactions based on CORUM complex membership
	def getPositiveAndNegativeInteractions(self):
		positive = set ([])
		negative = set([])
		prot2cluster = self.getProtToComplexMap()
		for protA in prot2cluster:
			for protB in prot2cluster:
				if protA == protB: continue
				edge = "\t".join(sorted([protA, protB]))
				if len(prot2cluster[protA] & prot2cluster[protB]) > 0:
					positive.add(edge)

				else:
					negative.add(edge)
		return positive, negative

	def overlap(self, a, b):
		tmpa = set(a)
		tmpb = set(b)
		overlap = math.pow((len(tmpa & tmpb)), 2) / (len(tmpa) * len(tmpb))
		return overlap

	# @author Florian Goebels
	# merges complexes which have an overlapp score > overlap_cutoff, and continues to merge until there is nothing left to merge
	def merge_complexes(self):
		merged = set()
		allComplexes = self.complexes.keys()
		newComplexes = {}

		for i in range(len(allComplexes)):
			compI = allComplexes[i]
			if compI in merged: continue
			candidates = [compI]
			toMerge = set([compI])
			while(len(candidates)>0):
				this_complex = candidates.pop()
				for j in range(i + 1, len(allComplexes)):
					compJ = allComplexes[j]
					if compJ in merged: continue
					if self.overlap(self.complexes[compI], self.complexes[compJ]) > self.overlap_cutoff:
						toMerge.add(compJ)
						candidates.append(compJ)
				merged.add(this_complex)


			if len(toMerge) > 0:
				(newName, newProts) = set(), set()
				for name in toMerge:
					newName.add(name)
					newProts.update(self.complexes[name])
				newComplexes[(",".join(map(str,newName)),)] = newProts
				merged.add(compI)
			else:
				newComplexes[compI] = self.complexes[compI]
				merged.add(compI)

		self.complexes = newComplexes

	def filter_complexes(self):
		todel = set()
		for comp in self.complexes:
			if len(self.complexes[comp]) < self.lb or len(self.complexes[comp]) > self.ub:
				todel.add(comp)
		for delComp in todel:
			del self.complexes[delComp]

	def remove_proteins(self, to_keep):
		for comp in self.complexes:
			self.complexes[comp] = self.complexes[comp] & to_keep
		self.filter_complexes()

	def getProtToComplexMap(self):
		out = {}
		for cluster in self.complexes:
			for prot in self.complexes[cluster]:
				if prot not in out: out[prot] = set([])
				out[prot].add(cluster)
		return out

	def getOverlapp(self, complexesB, cutoff = 0.5):
		out = 0
		for comp_ID_A in self.complexes.keys():
			protsA = self.complexes[comp_ID_A]
			matched = False
			for comp_ID_B in complexesB.complexes.keys():
				protsB = complexesB.complexes[comp_ID_B]
				if self.overlap(protsA, protsB) > cutoff:
					matched = True
					break

			if matched: out += 1
		return out

	def mmr(self, reference):
		matchingscores = {}
		for ref_complex in reference.complexes:
			if ref_complex not in matchingscores: matchingscores[ref_complex] = 0
			for complex in self.complexes:
				matchingscores[ref_complex] = max(matchingscores[ref_complex], self.overlap(self.complexes[complex],  reference.complexes[ref_complex]))

		mmr = matchingscores.values()
		mmr = sum(mmr)/len(mmr)
		return mmr

	def frac_match_comp(self, reference):
		out_simcoe = set([])
		out_overlap = set([])
		out_comb = set([])
		def simco(a,b):
			tmpa = set(a)
			tmpb = set(b)
			return len(tmpa&tmpb)/(min(len(tmpa), len(tmpb)))

		for complex in self.complexes:
			for ref_complex in reference.complexes:
				overlap_score = self.overlap(self.complexes[complex],  reference.complexes[ref_complex])
				simco_score = simco(self.complexes[complex],  reference.complexes[ref_complex])
				if overlap_score >= 0.25:
					out_overlap.add(complex)
				if simco_score > 0.5:
					out_simcoe.add(complex)
				if (simco_score+overlap_score)/2 >= 0.375:
					out_comb.add(complex)

		out_overlap = len(out_overlap)/len(self.complexes)
		out_simcoe = len(out_simcoe)/len(self.complexes)
		out_comb = len(out_comb) / len(self.complexes)
		return "%f\t%f\t%f" % (out_overlap, out_simcoe, out_comb)

	def sensitivity(self, reference):
		max_overlap_per_predicted_clustes = {}
		sum_of_all_cluster_sizes = 0
		for predicted_cluster in self.complexes:
			if predicted_cluster not in max_overlap_per_predicted_clustes: max_overlap_per_predicted_clustes[predicted_cluster] = 0
			sum_of_all_cluster_sizes += len(self.complexes[predicted_cluster])
			for reference_cluster in reference.complexes:
				overlap = len(self.complexes[predicted_cluster] & reference.complexes[reference_cluster])
				max_overlap_per_predicted_clustes[predicted_cluster] = max( max_overlap_per_predicted_clustes[predicted_cluster], overlap)
		max_overlap_per_predicted_clustes = sum(max_overlap_per_predicted_clustes.values())
		return max_overlap_per_predicted_clustes/sum_of_all_cluster_sizes


	def ppv(self, reference):
		n = len(self.complexes)
		m = len(reference.complexes)
		overlap_mat = np.zeros((n,m))
		for n, predicted_cluster in enumerate(self.complexes):
			for m, reference_cluster in enumerate(reference.complexes):
				overlap = len(self.complexes[predicted_cluster] & reference.complexes[reference_cluster])
				overlap_mat[n,m] = overlap
		return np.sum(overlap_mat.max(axis=0))/np.sum(overlap_mat)

	def acc(self, reference, sn = None, ppv = None):
		if sn == None: sn = self.sensitivity(reference)
		if ppv == None: ppv = self.ppv(reference)
		return math.sqrt(sn*ppv)

	def clus_sep(self, reference):
		n = len(self.complexes)
		m = len(reference.complexes)
		row_F = np.zeros((n, m))
		col_F = np.zeros((n, m))

		for i, compA in enumerate(self.complexes):
			for j, compB in enumerate(reference.complexes):
				overlap = len(self.complexes[compA] & reference.complexes[compB])
				row_F[i,j] = overlap
				col_F[i,j] = overlap
		print np.sum(row_F)

		row_F = np.nan_to_num(row_F/np.sum(row_F, axis=1, keepdims=True))
		col_F = np.nan_to_num(col_F/np.sum(col_F, axis=0, keepdims=True))
		sep = np.nan_to_num(np.sum(row_F*col_F))
		sep_co = sep/m
		sep_cl = sep/n
		return math.sqrt(sep_co*sep_cl)

	def clus_eval(self, ref):
		mmr = self.mmr(ref)
		ppv = self.ppv(ref)
		sn = self.sensitivity(ref)
		acc = self.acc(ref, sn, ppv)
		sep = self.clus_sep(ref)
		prc = self.frac_match_comp(ref)
		return mmr, prc, sn, ppv, acc, sep


# @author Florian Goebels
# Wrapper class for retrieving go annotation for a given taxid from the QuickGo webservice : https://www.ebi.ac.uk/QuickGO/
class QuickGO():
	# @author Florian Goebels
	# object constructor makes internet connection and downloads go annotations into the wd/go_files folder as taxid.go file
	# @param
	#		taxid of species that go annotation should be downloaded
	def __init__(self, taxid):
		self.taxid = taxid
		self.complexes = Clusters(need_to_be_mapped=False)
		self.get_GO_complexes()


	# @author Florian Goebels
	# reads in go flat file if gaf 20 format as protein to go annotation mapping (as dictonary)
	# @param
	#		taxid species for which go annotation should be read into memory
	def get_GO_complexes(self):
		go_to_prot_map = {}
		prot_to_go_map = {}
		quickgoURL = "http://www.ebi.ac.uk/QuickGO/GAnnotation?goid=GO:0043234&tax=%s&format=tsv&limit=1000000000&evidence=IDA,IPI,EXP" % (self.taxid)
		quickgoURL_FH = urllib2.urlopen(quickgoURL)
		quickgoURL_FH.readline()
		for line in quickgoURL_FH:
			line = line.rstrip()
			linesplit = line.split("\t")
			prot = linesplit[1]
			go_complex = linesplit[6]
			# Adding prot to go map
			if prot not in prot_to_go_map: prot_to_go_map[prot] = set([])
			prot_to_go_map[prot].add(go_complex)
			# Adding go to prot map
			if go_complex not in go_to_prot_map: go_to_prot_map[go_complex] = set([])
			go_to_prot_map[go_complex].add(prot)

		quickgoURL_FH.close()
		i = 0
		for go_complex in go_to_prot_map:
			self.complexes.addComplex(i, go_to_prot_map[go_complex])
			i+=1

	def get_complexes(self):
		return self.complexes



# @author Florian Goebels
# Class for handling and getting data from INPARANOID http://inparanoid.sbc.su.se/cgi-bin/index.cgi
class Inparanoid():
	# @author Florian Goebels
	# init function fetches inparanoid data and read it in for creating CORUM reference data
	# thus it always parses and download human to target species mapping
	# @param
	#		taxid taxid of given target species from which should be mapped to human
	#		inparanoid_cutoff cut off used for accepting ortholog mappings TODO include bootstrapping score as well
	#		foundProts list of proteins that were found in the MS experiments. If given, in cases of uncertain ortholog mapping (e.g. more than one mapping if inparanoid score == 1), ortholog mappings to found prots are prefered
	def __init__(self, taxid, inparanoid_cutoff=1, foundProts = set([])):
		self.taxid2Name = self.getTaxId2Namemapping()
		self.inparanoid_cutoff = inparanoid_cutoff
		self.foundProts = foundProts
		if taxid in self.taxid2Name:
			self.species = self.taxid2Name[taxid]
			xmldoc = self.getXML()
			self.orthmap, self.orthgroups = self.parseXML(xmldoc)
		else:
			print "Taxid:%s not supported" % taxid

	# @author Florian Goebels
	# mappes protein interactions to their respectiv ortholog counterpart
	# here it maps human reference protein interactions from CORUM to given target species
	# @param
	def mapEdges(self, edges):
		mapped_edges = set([])
		for edge in edges:
			protA, protB = edge.split("\t")
			if protA not in self.orthmap or protB not in self.orthmap: continue
			edge = "\t".join(sorted([self.orthmap[protA], self.orthmap[protB]]))
			mapped_edges.add(tuple(edge))
		return mapped_edges

	def mapComplexes(self, clusters):
		for clust in clusters.complexes:
			mapped_members = set([])
			for prot in clusters.complexes[clust]:
				if prot in self.orthmap:
					mapped_members.add(self.orthmap[prot])
#				else:
#					print "No map for %s" % prot
			clusters.complexes[clust] = mapped_members

	# @author Florian Goebels
	# get taxid to inparanoid name mapping from the inparanoid website
	def getTaxId2Namemapping(self):
		#TODO do not hard code
		url_str = "http://inparanoid.sbc.su.se/download/current/sequences/species.mapping.inparanoid8"
		url_FH = urllib2.urlopen(url_str)
		taxid2Name = {}
		for line in url_FH:
			line = line.rstrip()
			(taxid, name) = line.split(".fasta\t")
			taxid2Name[taxid] = name
		return taxid2Name

	# @author Florian Goebels
	# reads in ortholog mapping from the given inparanoid XML object
	# @param:
	#		xmldoc parsed Inparanoid xml doc
	def parseXML(self, xmldoc):
		protID2prot = {}
		orthgroups = []
		targetGenes = set([])
		for species in xmldoc.getElementsByTagName('species'): 
			for protlist in species.getElementsByTagName('genes'):
				for prot in protlist.getElementsByTagName('gene'):
					if species.getAttribute('NCBITaxId') != "9606": targetGenes.add(str(prot.getAttribute('protId')))
					protID2prot[prot.getAttribute('id')] = prot.getAttribute('protId')
		for orthgroup in xmldoc.getElementsByTagName('orthologGroup'):
			protsInGroup = set([])
			for protRef in orthgroup.getElementsByTagName('geneRef'):
				prot = str(protID2prot[protRef.getAttribute('id')])
				inparanoidScore = float(protRef.getElementsByTagName('score')[0].getAttribute('value'))
	#			TODO integrate bootstraping score
	#			bootstrapScore = protRef.getElementsByTagName('score')[1].getAttribute('value')
				if inparanoidScore >= self.inparanoid_cutoff:
					protsInGroup.add(prot)
			if len(protsInGroup)>1:
				orthgroups.append(protsInGroup)

		if len(self.foundProts) != 0:
			toDel = targetGenes - self.foundProts
			for i in range(len(orthgroups)):
				orthgroups[i] = orthgroups[i] - toDel
		outmap = {}
		outgroups = []
		for orthgroup in orthgroups:
			if len(orthgroup) == 2:
				protA, protB = orthgroup
				if protA in targetGenes:
					outmap[protB] = protA
				else:
					outmap[protA] = protB
				outgroups.append(orthgroup)
		return outmap, outgroups

	# @author Florian Goebels
	# fetches human to target Imparanoid xml doc from the inparanoid web server
	def getXML(self):
		url_str = 'http://inparanoid.sbc.su.se/download/current/Orthologs_OrthoXML/%s/%s-%s.orthoXML'
		first = self.species
		second = "H.sapiens"	
		if self.species > "H.sapiens":
			first = "H.sapiens"
			second = self.species
		url_str = url_str % (first, first, second)
		xml_str = urllib2.urlopen(url_str).read()
		xmldoc = minidom.parseString(xml_str)
		return xmldoc

	def readTable(self, tableF):
		tableFH = open(tableF)
		targetGenes = set([])
		orthgroups = []
		def getids(ids_raw):
			out = []
			for i in range(int(len(ids_raw)/2)):
				this_id = ids_raw[i*2]
				this_score = float(ids_raw[(i*2)+1])

				if this_score < self.inparanoid_cutoff:continue
				out.append(this_id)
			return out
		tableFH.readline()
		for line in tableFH:
			orthID, score, target_orth, source_orth = line.rstrip().split("\t")
			targetids = getids(target_orth.split())
			sourceids = getids(source_orth.split())
			group = []
			group.extend(targetids)
			group.extend(sourceids)
			orthgroups.append(set(group))
			targetGenes |= set(targetids)
		tableFH.close()

		if len(self.foundProts) != 0:
			toDel = targetGenes - self.foundProts
			for i in range(len(orthgroups)):
				orthgroups[i] = orthgroups[i] - toDel
		outmap = {}
		outgroups = []
		for orthgroup in orthgroups:
			if len(orthgroup) == 2:
				protA, protB = orthgroup
				if protA in targetGenes:
					outmap[protB] = protA
				else:
					outmap[protA] = protB
				outgroups.append(orthgroup)
		self.orthmap, self.orthgroups =  outmap, outgroups