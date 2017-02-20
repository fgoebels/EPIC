#import glob
#import sys
#import urllib
#import urllib2
#import os
#import re

#taxoID = '6239'

#taxoIDspeciesDic = {'3702':'Arabidopsis_thaliana', '6239':'Caenorhabditis_elegans', '7955':'Danio_rerio', 
                                    #'7227':'Drosophila_melanogaster','562':'Escherichia_coli','9606':'Homo_sapiens',
                                    #'10090':'Mus_musculus','10116':'Rattus_norvegicus','4932':'Saccharomyces_cerevisiae'} 

#taxoIDurl = {'6239':'http://www.uniprot.org/uniprot/?query=taxonomy%3A6239&sort=score&columns=id,genes(ORF)&format=tab'}

#Caenorhabditis_elegans
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A6239&sort=score&columns=id,genes(ORF)&format=tab' #ok

#Arabidopsis_thaliana
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A3702&sort=score&columns=id,genes(OLN)&format=tab' #ok

#Danio_rerio
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A7955&sort=score&columns=id,database(Ensembl)&format=tab' #ok

#Drosophila_melanogaster
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A7227&sort=score&columns=id,database(FlyBase)&format=tab' #ok

#Escherichia_coli
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A562&sort=score&columns=id,database(Ensembl)&format=tab' 

#Homo_sapiens
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A9606&sort=score&columns=id,database(Ensembl)&format=tab'

#Mus_musculus
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A10090&sort=score&columns=id,database(Ensembl)&format=tab'

#Rattus_norvegicus
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A10116&sort=score&columns=id,database(Ensembl)&format=tab'

#Saccharomyces_cerevisiae
#url = 'http://www.uniprot.org/uniprot/?query=taxonomy%3A4932&sort=score&columns=id,genes&format=tab' #ok

#response = urllib2.urlopen(taxoIDurl[taxoID])

#html = response.readlines() #return everything in a list, each item in the list is a line in original file

#unipront_geneNames_dic = {}

#for item in html[1:]:
    
    #items = item.split("\t")
    #uniprot = items[0]
    
    #geneNames = items[1].strip()
    
    #genes_list = re.split('[;\s\|]', geneNames)
    
    #new_list = list(genes_list) #make a new list to clone original list, otherwise, the code wont work.
    
    #for gene_names in genes_list:
	
	#if "CELE" in gene_names:
	    #new_list.remove(gene_names)
    
    #new_list = filter(None, new_list) #remove the empty item from the list
	    
    #if len(new_list) > 0 :
	#unipront_geneNames_dic[uniprot] = new_list
    
#for key, value in unipront_geneNames_dic.iteritems():
    #print key, value
    #print len(value)

from __future__ import division
import CalculateCoElutionScores as CS


def main():
    
    testing = CS.Genemania('6239')
    testing.map_proteinNames()










if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		pass