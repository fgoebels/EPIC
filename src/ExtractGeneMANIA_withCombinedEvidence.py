import glob
import sys
from collections import defaultdict


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

# all functional evidence codes in GeneMANIA, excluding "Physical" and "complexes" to eliminate circularity
functionalEvidences = ['Co-expression', 'Genetic', 'Other', 'Predicted', 'Shared']

# read files of GeneMANIA database and add features to 

# create the 2-levels' master dictionary,
# first level key is evidence name
# second level key is sorted protein-pair
evidenceDictionary = {} 

for f_evidence in functionalEvidences:
    secondaryEvidenceDic = defaultdict(list)
            
    for filename in glob.glob('*.txt'):
        
        if filename.startswith(f_evidence):
                        
            print "Processing: %s" % (filename)
            fh = open(filename)
        
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
    geneMANIAOutput.write(edge)
    
    for i in range (0,len(proteinsPairDictionary[edge])):
        geneMANIAOutput.write("\t")
        geneMANIAOutput.write(str('%.5f' % proteinsPairDictionary[edge][i]))
    
    geneMANIAOutput.write("\n")
        
geneMANIAOutput.close()