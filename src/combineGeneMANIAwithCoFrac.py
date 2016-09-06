import sys

# @ author Lucas Ming Hu
# Combine GeneMANIA functional evidence with 
# co-fractionation experimental data.
# commanda line: python combineGeneMANIAwithCoFrac.py GeneMANIAFile coFractionationFile


proteinsPairGeneMANIADictionary = {}
noItems = 0;

noItemsFile1 = 0;        
with open(sys.argv[1]) as myfile:    
    for line in myfile:
        line = line.strip()
        items = line.split("\t")
        edge = "\t".join(sorted([items[0], items[1]]))
        noItemsFile1 = len(items[2:]) 
        proteinsPairGeneMANIADictionary[edge] = '\t'.join(items[2:])
myfile.close()
        
noItemsFile2 = 0
with open(sys.argv[2]) as myfile:
    for line in myfile:
        line = line.strip()
        items = line.split("\t")
        edge = "\t".join(sorted([items[0], items[1]]))
        value = '\t'.join(items[2:])
        noItemsFile2 = len(items[2:]) 
        if edge in proteinsPairGeneMANIADictionary:
            proteinsPairGeneMANIADictionary[edge] = proteinsPairGeneMANIADictionary[edge] + "\t" + value
myfile.close()

noItemsInBothFiles = noItemsFile1 + noItemsFile2
                
geneMANIAOutput = open((sys.argv[1])[:-4] + "_experimental_data_GeneMANIAFeatures.txt", "w")
for edge in proteinsPairGeneMANIADictionary:
    items = proteinsPairGeneMANIADictionary[edge].split("\t")
    if len(items) == noItemsInBothFiles:
        geneMANIAOutput.write(edge + "\t" + proteinsPairGeneMANIADictionary[edge] + "\n")
geneMANIAOutput.close()

            