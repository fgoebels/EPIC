#!/usr/bin/env python
from __future__ import division
import numpy as np
import random
import utils 
import sys
import urllib
from xml.dom import minidom
import os
import wget
import CalculateCoElutionScores as CalcS


# @author Florian Goebels
# place holder function for debugging
def main():
	(targetSpecies, elutionFile, outF) = sys.argv[1:]

	elutionData = CalcS.ElutionData(elutionFile)
	refProts = set(elutionData.prot2Index.keys())

	inparanoid = ""
	if targetSpecies != "9606":
		inparanoid = Inparanoid(targetSpecies, foundProts = refProts)
	gos = QuickGO(targetSpecies)
	corum = CORUM()
	reference = Goldstandard_from_CORUM(corum, inparanoid, ratio=1, found_prots = refProts)
	print len(reference.goldstandard_positive)
	print len(reference.goldstandard_negative)
	print reference.goldstandard_positive
#	print reference.goldstandard_negative

# @author Florian Goebels
# Class for creating reference set from CORUM
class Goldstandard_from_CORUM():
	# @author Florian Goebels
	# Constructor for creating and calculating refrence data set from CORUM complex data
	# @param
	#		complexes CORUM complex object Human complexes
	#		orthmap Inparanoid object for mapping Human CORUM complexes to target species
	#		ratio ratio for negative to postivie with |negative| = ratio * |positive|
	#		found_prots proteins identified via MS cofractionation either as list
	def __init__(self, complexes, orthmap="", ratio = 5, found_prots = ""):
		self.complexes = complexes
		self.ratio = ratio
		self.found_prots = found_prots
		self.goldstandard_positive = set([])
		self.goldstandard_negative = set([])
		self.getPositiveAndNegativeInteractions()
		if orthmap != "":
			self.orthmap = orthmap
			self.mapReferenceData()
		self.makeReferenceDataSet()

	# @author Florian Goebels
	# creats all possible positive and negative protein interactions based on CORUM complex membership
	def getPositiveAndNegativeInteractions(self):
		prot2cluster = self.complexes.getProtToComplexMap()
		for protA in prot2cluster:
			for protB in prot2cluster:
				if protA == protB: continue
				edge = tuple(sorted([protA, protB]))
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
	# creates reference data set with lable information and set ration of positive and negative protein interactions
	def makeReferenceDataSet(self):
		def remove_not_found_prots(protlist, foundprots):
			out = set([])
			for (protA, protB) in protlist:
				if protA not in foundprots or protB not in foundprots: continue
				out.add((protA, protB))
			return out
		self.goldstandard_positive = remove_not_found_prots(self.goldstandard_positive, self.found_prots)
		self.goldstandard_negative = remove_not_found_prots(self.goldstandard_negative, self.found_prots)
		reference = set([])
		pos = set([])
		neg = set([])
		for (protA, protB) in self.goldstandard_positive:
			pos.add((protA, protB, "positive"))

		for (protA, protB) in random.sample(self.goldstandard_negative, len(self.goldstandard_positive)*self.ratio):
			neg.add((protA, protB, "negative"))

		reference = pos
		reference.update(neg)
		self.goldstandard = reference


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
	def __init__(self, lb = 3, ub=30, overlap_cutoff=0.5, srouce_species = "Human"):
		self.corum_file = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2]) + os.sep + "data" + os.sep +"corum.txt"
		self.srouce_species = srouce_species
		self.ub = ub
		self.lb = lb
		self.overlap_cutoff = overlap_cutoff
		self.getCORUM()
		self.complexes = {}
		self.readCORUM()
		self.filterComplexes()
		self.mergeComplexes()

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
			corum_url = "http://mips.helmholtz-muenchen.de/genre/proj/corum/allComplexesCore.csv"
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
		complexes2prot = utils.readData(self.corum_file, np.array([0]), np.array([3,4]), sep=";")
		for comp in complexes2prot:
			for anno in complexes2prot[comp]:
				(species, prots) = anno
				if species != self.srouce_species: continue
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
			overlap = len(tmpa & tmpb)/len(tmpa | tmpb)
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
		for (protA, protB) in edges:
			if protA not in self.orthmap or protB not in self.orthmap: continue
			mapped_edges.add(tuple(sorted([self.orthmap[protA], self.orthmap[protB]])))
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
				inparanoidScore = protRef.getElementsByTagName('score')[0].getAttribute('value')
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
		

if __name__ == "__main__":
        try:
                main()
        except KeyboardInterrupt:
                pass
