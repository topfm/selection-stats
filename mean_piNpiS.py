from collections import defaultdict
import sys
import numpy

##check for correct args
if len(sys.argv) != 2:
    print("Usage: mean_piNpiS.py selectionStats_piNpiS.txt")
    sys.exit(0)

statsfile = sys.argv[1]


##read genes into file with dict = {gene1:[val1,val2,val3], gene2:[val1,val2,val3]}
##Here, I am skipping any piNpiS value of "None", which indicates that the codon contains an "N" (instead of ACTG)
dd = defaultdict(list)
avDict = {}
with open(statsfile, "r") as infile:
    next(infile)
    for line in infile:
        if "None" in line:
            pass
        else:
            info = line.strip().split("\t")
            gene= info[0]
            gene = gene.split("_")
            piNpiS = info[-2]
            if len(gene) > 3:
                gene = gene[0] + "_" + gene[1]
            else:
                gene = gene[0]
            dd[gene].append(piNpiS)
        

##Averaging each set of values for each key in the dictionary
##This will give an average piNpiS value for each gene!
for k, v in dd.items():
    v = [(float(i)) for i in v]
    mean = sum(v) / len(v)
    avDict[k] = round(mean, 2)


##Writing the output file
outfile = statsfile.replace(".txt", "_average.csv") 

with open(outfile, "w") as f:
     for key in avDict.keys():
        f.write("%s, %s\n" %(key, avDict[key]))
   




