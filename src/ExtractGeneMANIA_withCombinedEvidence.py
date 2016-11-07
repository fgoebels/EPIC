import glob
import sys
from collections import defaultdict
from bs4 import BeautifulSoup
import urllib
import urllib2
import os 

def catchFile(species): 
    speciesList = ['Danio_rerio']
    if species not in speciesList:
        return None 

    urlbase = 'http://genemania.org/data/current'
    speciesURL = os.path.join(urlbase, species)
    r = urllib.urlopen(speciesURL).read()
    soup = BeautifulSoup(r)

    table = soup.find('table')

    allcell = []
    for row in table.find_all('tr'):
        for col in row.find_all('td'):
            allcell.append(col.getText())

    #filtering 
    result = [] 
    for c in allcell:
        if '.txt' in c:
            result.append(os.path.join(speciesURL,c))

    return result



# @ author Lucas Ming Hu
# Extract functional evidence from GeneMANIA database
# Combine evidences from the same source into a single line

#create a dictionary to store all input protein-protein pairs.
proteinsPairDictionary = defaultdict(list)

with open(sys.argv[1]) as myfile:    
    for line in myfile:
        items = line.split()
        edge = "\t".join(sorted([items[0], items[1]]))
        proteinsPairDictionary[edge] = []

# @author: Lucas Ming Hu        
# a helper function to get the average of the GeneMANIA scores 
# for each line of evidence
def average(secondaryEvidenceDic):
    resultDict = defaultdict(float)
    
    for key in secondaryEvidenceDic:
        resultDict[key] = sum(secondaryEvidenceDic[key]) * 1.0 / len(secondaryEvidenceDic[key])
    
    return resultDict 

# all functional evidence codes in GeneMANIA, excluding "Physical" and "complexes" and "Predicted" to eliminate circularity
functionalEvidences = ['Co-expression', 'Genetic', 'Other', 'Shared']

# read files of GeneMANIA database and add features to 

# read online database - by species name
fps = catchFile('Danio_rerio')

# create the 2-levels' master dictionary,
# first level key is evidence name
# second level key is sorted protein-pair
evidenceDictionary = {} 

for f_evidence in functionalEvidences:
    secondaryEvidenceDic = defaultdict(list)
            
    for fp in fps:
        filename = fp.split('/')[-1]

        if filename.startswith(f_evidence):
                        
            print "Processing: %s" % (filename)
            fh = urllib2.urlopen(fp) 
            
            for line in fh: 
                proteinA, proteinB, score = line.split()
                edge = "\t".join(sorted([proteinA, proteinB]))

                if edge in proteinsPairDictionary: 
                    secondaryEvidenceDic[edge].append(float(score))
                    
            fh.close()
    
    evidenceDictionary[f_evidence] = average(secondaryEvidenceDic)      


# combine all information from the 2-level dictionary "evidenceDictionary" into "proteinsPairDictionary"
# the key for "proteinsPairDictionary" is protein-pair, and the each line is a combined value

for key1 in evidenceDictionary:
    for key2 in proteinsPairDictionary:
        if key2 in evidenceDictionary[key1]:
            proteinsPairDictionary[key2].append(evidenceDictionary[key1][key2])   
        else:       
            proteinsPairDictionary[key2].append(0)
                        

geneMANIAOutput = open((sys.argv[1])[:-4] + "_GeneMANIAFeatures_combinedEvidence.txt", "w")
                   
for edge in proteinsPairDictionary:
    
    if sum(proteinsPairDictionary[edge]) > 0:
        geneMANIAOutput.write(edge)
            
        for i in range (0,len(proteinsPairDictionary[edge])):
            geneMANIAOutput.write("\t")
            geneMANIAOutput.write(str('%.5f' % proteinsPairDictionary[edge][i]))
    
        geneMANIAOutput.write("\n")
        
geneMANIAOutput.close()