import glob
import sys


# @ author Lucas Ming Hu
# Extract functional evidence from GeneMANIA database  

#create a dictionary to store all input protein-protein pairs.
proteinsPairDictionary = {}

with open(sys.argv[1]) as myfile:    
    for line in myfile:
        proteinA, proteinB= line.split()
        edge = "\t".join(sorted([proteinA, proteinB]))
        proteinsPairDictionary[edge] = ""

# @author: Lucas Ming Hu
# helper function for checking if the dictionary of protein-protein pairs
# has the right number of scores after reading each file, the missing value
# is filled by zero.
def fillZero(proteinsPairDictionary, proteinsInFile):
    proteins_I_have = set(proteinsPairDictionary.keys())
    proteins_wih_no_evidence = proteins_I_have - proteinsInFile
    for protein in proteins_wih_no_evidence:
        proteinsPairDictionary[protein] += "\t0"

#all functional evidence codes in GeneMANIA
functionalEvidences = ['Co-expression', 'Genetic', 'other', 'Physical', 'Predicted', 'Shared']

#read files of GeneMANIA database and add features to 
for filename in glob.glob('*.txt'):
    has_evidence = False
    for f_eveidence in functionalEvidences:
        if f_eveidence in filename:
            has_evidence = True
    if not has_evidence: continue
    print "Processing: %s" % (filename)
    fh = open(filename)
    edges_in_file = set([])
    for line in fh:  
        proteinA, proteinB, score = line.split()
        edge = "\t".join(sorted([proteinA, proteinB]))
        edges_in_file.add(edge)
        if edge not in proteinsPairDictionary: continue
        proteinsPairDictionary[edge] += "\t%s" % score
    fh.close()
    fillZero(proteinsPairDictionary, edges_in_file)
    edges_in_file = set([])    


geneMANIAOutput = open((sys.argv[1])[:-4] + "_GeneMANIAFeatures.txt", "w")
                            
for edge in proteinsPairDictionary:
    geneMANIAOutput.write("%s%s\n" % (edge, proteinsPairDictionary[edge]))
geneMANIAOutput.close()