#! /usr/bin/env python

import sys
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

####
# - Takes a whole genome alignment file (RGA) and an annotation file (in bed format)
#   and splits the alignment into individual gene alignments.
# - Will name gene alignments by "Name" if present in info column of annotation file
#   otherwise will name them by "ID".
# - Gene alignments will be written to a directory called "indv_genes/"
#   ** Note: features that are not of the type "CDS" will be written into a directory called "indv_genes/non_CDS/" 
####

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: alignmentToGenes.py genome_alignment.fasta annotation.bed")
    sys.exit(0)

alnfile = sys.argv[1]
annotfile = sys.argv[2]

# make output directory
subprocess.call(["mkdir","indv_genes"])
subprocess.call(["mkdir","indv_genes/non_CDS"])

# read annotations into a dictionary with structure {"gene":[start,stop]}
genes = {}
with open(annotfile,"r") as f:
    for line in f:
        info = line.strip().split("\t")
        start = int(info[1])
        stop = int(info[2])
        typ = info[7]
        annot = info[9].split(";")
        if "Name" not in info[9]:
            for x in annot:
                if x.startswith("ID"):
                    name = x.split("=")[1]
        else:
            for x in annot:
                if x.startswith("Name"):
                    name = x.split("=")[1]
        genes[name] = [start,stop,typ]

# parse genome alignment, for each sequence read through genes dictionary,
# create new Seq object with gene sequence and append gene sequence to new dictionary
# with structure {"gene":[Seq1,Seq2]}
indvalns = defaultdict(list)
for genome in SeqIO.parse(alnfile,"fasta"):
    for gene,coord in genes.items():
        alnname = genome.id+"_"+gene
        geneseq = SeqRecord(Seq(str(genome.seq[coord[0]:coord[1]])),\
                id=alnname,\
                name=alnname,\
                description=alnname)
        indvalns[gene].append(geneseq)

# write individual alignments
for gene,alns in indvalns.items():
    if genes[gene][2] == "CDS":
        newfile = open("indv_genes/"+gene+".fasta","w")
        SeqIO.write(alns,newfile,"fasta")
        newfile.close()
    else:
        newfile = open("indv_genes/non_CDS/"+gene+".fasta","w")
        SeqIO.write(alns,newfile,"fasta")
        newfile.close()
