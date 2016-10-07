#!/usr/bin/python
from __future__ import division
import numpy as np
import random
import sys
import urllib
from xml.dom import minidom
import os
import wget
import re


# @author Florian Goebels
# place holder function for debugging
def main():
	(targetSpecies, elutionFile, outF) = sys.argv[1:]

	elutionData = CalcS.ElutionData(elutionFile)
	refProts = set(elutionData.prot2Index.keys())

	reference = Goldstandard_from_CORUM("9606", found_prots = refProts)
#	print len(reference.goldstandard_positive)
#	print len(reference.goldstandard_negative)
#	print reference.goldstandard_positive
#	print reference.goldstandard_negative

class GS_from_PPIs():
	def __init__(self, refF, found_ppis=""):
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
				if edge not in found_ppis: continue
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
		self.goldstandard_positive, self.goldstandard_negative = self.readGS(gsF)

	def readGS(self, gsF):
		negative = set([])
		positive = set([])
		gsFH = open(gsF)
		prot2clusters = {}
		i = 0
		for line in gsFH:
			line = line.rstrip()
			prots = line.split("\t")
			for prot in prots:
					if prot not in prot2clusters: prot2clusters[prot] = set([])
					prot2clusters[prot].add(i)
			i+=1
		prots = prot2clusters.keys()
		for i in range(len(prots)):
			protA = prots[i]
			clustA = prot2clusters[protA]
			for j in range(i+1, len(prots)):
				protB = prots[j]
				clustB = prot2clusters[protB]
				edge = "\t".join(sorted([protA, protB]))
				if len(clustA & clustB) > 0:
					positive.add(edge)
				else:
					negative.add(edge)
		return positive, negative

# @author Florian Goebels
# Class for creating reference set from CORUM
class Goldstandard_from_CORUM():
	# @author Florian Goebels
	# Constructor for creating and calculating refrence data set from CORUM complex data
	# @param
	#		complexes CORUM complex object Human complexes
	#		orthmap Inparanoid object for mapping Human CORUM complexes to target species
	#		found_prots proteins identified via MS cofractionation either as list
	def __init__(self, targetSpecies, found_prots = "", source_species_regex = "(Human|Mammalia)"):
		inparanoid = ""
		if targetSpecies != "9606" and "Human" in source_species_regex:
			inparanoid = Inparanoid(targetSpecies, foundProts=found_prots)
		corum = CORUM(source_species_regex = source_species_regex)
		self.complexes = corum
		self.found_prots = found_prots
		self.goldstandard_positive = set([])
		self.goldstandard_negative = set([])
		self.getPositiveAndNegativeInteractions()
		if inparanoid != "":
			self.orthmap = inparanoid
			self.mapReferenceData()

	# @author Florian Goebels
	# creats all possible positive and negative protein interactions based on CORUM complex membership
	def getPositiveAndNegativeInteractions(self):
		prot2cluster = self.complexes.getProtToComplexMap()
		for protA in prot2cluster:
			for protB in prot2cluster:
				if protA == protB: continue
				edge = "\t".join(sorted([protA, protB]))
				if len(prot2cluster[protA] & prot2cluster[protB]) > 0:
					self.goldstandard_positive.add(edge)

				else:
					self.goldstandard_negative.add(edge)
		

	# @author Florian Goebels
	# maps reference protein complexes from CORUM source species (Human) to target species
	def mapReferenceData(self):
		self.goldstandard_positive = self.orthmap.mapEdges(self.goldstandard_positive)
		self.goldstandard_negative = self.orthmap.mapEdges(self.goldstandard_negative)

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
	def __init__(self, lb = 1, ub=50, overlap_cutoff=0.5, source_species_regex = "(Human|Mammalia)"):

		# static regex for identifying valid bochemical evidences codes
		self.biochemical_evidences_regex ="MI:(2193|2192|2191|2197|2195|2194|2199|2198|0807|0401|0400|0406|0405|0404|0089|0084|0081|0007|0006|0004|0513|1029|0979|0009|0008|0841|1312|2188|2189|0411|0412|0413|0928|0415|0417|0098|0729|0920|0921|0603|0602|0605|0604|0402|0095|0096|0606|0091|0092|1142|1145|1147|0019|1309|0696|0697|0695|0858|0698|0699|0425|0424|0420|0423|0991|0990|0993|0992|0995|0994|0997|0996|0999|0998|1028|1011|1010|1314|0027|1313|0029|0028|0227|0226|0225|0900|0901|0430|0434|0435|1008|1009|0989|1004|1005|0984|1007|1000|0983|1002|1229|1087|1325|0034|0030|0031|0972|0879|0870|1036|0678|1031|1035|1034|0676|0440|1138|1236|0049|0048|1232|0047|1137|0419|0963|1026|1003|1022|0808|0515|0514|1187|0516|0511|1183|0512|0887|0880|0889|0115|1006|1249|0982|0953|1001|0508|0509|0657|0814|1190|1191|0813|0066|0892|0899|1211|0108|1218|1352|1354|0949|0946|0947|0073|0071|1019|2168|0700|2167|1252|1017|0276|1189|1184)"
		self.corum_file = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2]) + os.sep + "data" + os.sep +"corum.txt"
		self.source_species_regex = source_species_regex
		self.ub = ub
		self.lb = lb
		self.overlap_cutoff = overlap_cutoff
		self.getCORUM()
		self.complexes = {}
		self.readCORUM()
		self.filterComplexes()
#		self.mergeComplexes()

	# @author Florian Goebels
	# creates protein to each complex asocciation mapping,
	# which is used as helper function to find proteins that share at least one corum complex.
	# This funcrion is used by Goldstandard_from_CORUM class to create the reference data set
	def getProtToComplexMap(self):
		out = {}
		for cluster in self.complexes:
			for prot in self.complexes[cluster]:
				if prot not in out: out[prot]=set([])
				out[prot].add(cluster)
		return out

	# @author Florian Goebels
	# downloads current version of corum and safe it to wd/data folder as corum.txt
	def getCORUM(self):
		if(not os.path.exists(self.corum_file)):
			corum_url = "http://mips.helmholtz-muenchen.de/genre/proj/corum/allComplexes.csv"
			wget.download(corum_url, self.corum_file)

	# @author Florian Goebels
	# filteres complexes with |size| > lb or |size| > ub
	# @param
	#		lb lower bound complex should have at least lb members
	#		ub upper bound complex should have at most  ub members
	def filterComplexes(self):
		todel = set()
		for comp in self.complexes:
			if len(self.complexes[comp]) < self.lb or len(self.complexes[comp]) > self.ub:
				todel.add(comp)
		for delComp in todel:
			del self.complexes[delComp]
			

	# @author Florian Goebels
	# reads in CORUM from flat file
	def readCORUM(self):
		complexes2prot = readData(self.corum_file, np.array([0]), np.array([3,4,6]), sep=";")
		for comp in complexes2prot:
			for anno in complexes2prot[comp]:
				(species, prots, evidence) = anno
				# bool(...) returns true if evidence code is found => not bool(...) is true if not valid evidence is found, and thus skip this complex
				if not bool(re.search(self.biochemical_evidences_regex, evidence)):
					continue
				if not bool(re.search(self.source_species_regex, species)): continue
				prots = prots.split(",")
				self.complexes[comp] = prots

	# @author Florian Goebels
	# merges complexes which have an overlapp score > overlap_cutoff, and continues to merge until there is nothing left to merge
	def mergeComplexes(self):
		merged = set()
		allComplexes = self.complexes.keys()
		newComplexes = {}
		def overlap(a,b):
			tmpa = set(a)
			tmpb = set(b)
			overlap = len(tmpa & tmpb)/min(len(tmpa),len(tmpb))
			return overlap

		for i in range(len(allComplexes)):
			toMerge = set()
			compI = allComplexes[i]
			if compI in merged: continue
			for j in range(i+1, len(allComplexes)):
				compJ = allComplexes[j]
				if compJ in merged: continue
				if overlap(self.complexes[compI],self.complexes[compJ]) > self.overlap_cutoff:
					toMerge.add(compJ)
					merged.add(compJ)
				if len(toMerge)>0:
					toMerge.add(compI)
					(newName,  newProts) = set(), set()
					for name in toMerge:
						newName.add(name[0])
						newProts.update(self.complexes[name])
					newComplexes[(",".join(newName),)] = newProts
					merged.add(compI)
				else:
					newComplexes[compI] = self.complexes[compI]
					merged.add(compI)

		self.allComplexes = newComplexes

# @author Florian Goebels
# Wrapper class for retrieving go annotation for a given taxid from the QuickGo webservice : https://www.ebi.ac.uk/QuickGO/
class QuickGO():
	# @author Florian Goebels
	# object constructor makes internet connection and downloads go annotations into the wd/go_files folder as taxid.go file
	# @param
	#		taxid of species that go annotation should be downloaded
	def __init__(self, taxid):
		self.godir = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2]) + os.sep + "go_files" + os.sep +"%s.go"
		self.goMap = {}
		self.loadGO(taxid)

	# @author Florian Goebels
	# return go annotaiton for a given prot identifier (here Uniprot)
	# @param
	#		prot protein for which go annotation should be returned
	def getGOforProt(self, prot):
		out = set([])
		if prot in self.goMap: out = self.goMap[prot]
		return out

	# @author Florian Goebels
	# checks if go file exists locally, if not download go from quick go service
	# @param
	#		taxid taxid for which go annotion should be downloaded
	def loadGO(self, taxid):
		if(not os.path.exists(self.godir % ( taxid))):
			self.getGO(taxid)
		self.readGO(taxid)

	# @author Florian Goebels
	# reads in go flat file if gaf 20 format as protein to go annotation mapping (as dictonary)
	# @param
	#		taxid species for which go annotation should be read into memory
	def readGO(self, taxid):
		goF = self.godir % (taxid)
		goFH = open(goF)
		for line in goFH:
			line = line.rstrip()
			linesplit = line.split("\t")
			if len(linesplit) < 5: continue
			protID = linesplit[1]
			goID = linesplit[4]
			if protID not in self.goMap: self.goMap[protID] = set([])
			self.goMap[protID].add(goID)
		goFH.close()

	# @author Florian Goebels
	# automaticallys downloads go anntotation for a given tax id
	# TODO create switch to select which level of go annotaiton to by downloaded IEA and Manula, right now it only gets mahual annotation
	# @param
	#		taxid species for which go annotation should be downloaded
	def getGO(self, taxid):
		quickgoURL = "http://www.ebi.ac.uk/QuickGO/GAnnotation?tax=%s&format=gaf&limit=1000000000&evidence=IMP,CIGI,CIPI,CIDA,CIEP,CEXP,CISS,CNAS,CTAS,CND,CIC,CRCA,CIBA,CIBD,CIKR,CIRD,CISA,CISM,CISO,CIGC"
		wget.download(quickgoURL % (taxid), self.godir % ( taxid))
		
		
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
		self.species = self.taxid2Name[taxid]
		self.inparanoid_cutoff = inparanoid_cutoff
		self.foundProts = foundProts

		xmldoc = self.getXML()
		self.orthmap, self.orthgroups = self.parseXML(xmldoc)

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

	# @author Florian Goebels
	# get taxid to inparanoid name mapping from the inparanoid website
	def getTaxId2Namemapping(self):
		#TODO do not hard code
		url_str = "http://inparanoid.sbc.su.se/download/current/sequences/species.mapping.inparanoid8"
		url_FH = urllib.urlopen(url_str)
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
		xml_str = urllib.urlopen(url_str).read()
		xmldoc = minidom.parseString(xml_str)
		return xmldoc
		

def readData(dataF, keyColNumbers, valueRowNumbers, primKeyMap = "", header = True, sep="\t"):
        out = {}
        dataFH = open(dataF)
        if header: dataFH.readline()
        for line in dataFH:
                line = line.rstrip()
                lineSplit = np.array(line.split(sep))
#               if len(lineSplit)<2: continue
                key = tuple(lineSplit[keyColNumbers])
                primkey = tuple([key[0]])
                if primKeyMap != "" and primkey in primKeyMap:
                        if len(primKeyMap[primkey]) ==1:
                                newPrimKey = primKeyMap[primkey][0][0]
                                key = list(key[1:])
                                key.insert(0, newPrimKey)
                                key = tuple(key)
#                       else:
#                               print "no uniq map for " + primkey[0]
#               if primKeyMap != "" and primkey not in primKeyMap:
#                       print "no mapping found for " + primkey[0]
                value = tuple(lineSplit[valueRowNumbers])
                if key not in out: out[key] = []
                if value not in out[key]: out[key].append(value)
        dataFH.close()
        return out

if __name__ == "__main__":
		try:
				main()
		except KeyboardInterrupt:
				pass
