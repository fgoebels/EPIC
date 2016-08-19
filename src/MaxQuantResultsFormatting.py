import sys

# @ author Lucas Ming Hu
# returns correct format for EPIC from MaxQuant results file: proteinGroups.txt
# command line: python MaxQuantResultsFormatting.py proteinGroups.txt outputDir
# results files are: SpectralCountMS2EpicInput.txt and ItensityMS1EpicInput.txt
            
with open(sys.argv[1]) as myfile:
    
    content = myfile.readlines()[0:1]
    contentString = "".join(content)
    header = contentString.split("\t")
    
    for i in range(len(header) - 1):
        
        if ("Fasta headers") == header[i]:
            uniprotProteinNameIndex = i
            
        if ("Razor") in header[i]:
            MS2FractonsStartsIndex = i + 1
            
        if ("Sequence coverage [%]") == header[i]:
            MS2FractonsEndsIndex = i - 1
            
        if ("Intensity") == header[i]:
            MS1FractionStartsIndex = i + 1
            
        if ("MS/MS Count") == header[i]:
            MS1FractionEndsIndex = i - 1

spectralCount = open(sys.argv[2] + "SpectralCountMS2EpicInput.txt", "w")
intensity = open(sys.argv[2] + "ItensityMS1EpicInput.txt", "w")
    
with open(sys.argv[1]) as myfile:     
    for line in myfile:
        eachLine = line.split("\t")
        eachLineItem = eachLine[uniprotProteinNameIndex].split("|")
        
        if len(eachLineItem) > 1:    
            uniprotName = eachLineItem[1]
            intensity.write(uniprotName + "\t")
            spectralCount.write(uniprotName + "\t")
            
            for i in range(MS2FractonsStartsIndex, MS2FractonsEndsIndex):
                spectralCount.write(eachLine[i] + "\t")
            spectralCount.write(eachLine[MS2FractonsEndsIndex] + "\n")
            
            for i in range(MS1FractionStartsIndex, MS1FractionEndsIndex):
                intensity.write(eachLine[i] + "\t")
            intensity.write(eachLine[MS1FractionEndsIndex] + "\n")
            
