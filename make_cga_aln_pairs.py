#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter
import glob


if len(sys.argv) != 2:
	print("Usage: script.py <geneList.txt>")
	sys.exit(0)

def get_ids(idFile):
	idList = []
	with open(idFile,"r") as f:
		for line in f:
			idList.append(line.strip())
	return idList

###################################################
# iterate through each sequence from each gene aln
# 	for each seq in the aln:
# 		pair w each other seq (no repeats)
#		write fasta (with two seq)
###################################################

def make_alns(genes):
	for index,line in enumerate(genes):
		gene = str(line.strip())
		pairs = []
		alnName = "./core_gene_alns/"+gene+".fasta"
		allseq = []
		for seq_record in SeqIO.parse(alnName, "fasta"):
			allseq.append(seq_record)   
		for seq1 in allseq:
			seq1name = str(seq1.id.split("_")[0])
			for seq2 in allseq:
				seq2name = str(seq2.id.split("_")[0])
				if seq1name != seq2name:
					pairseqname = seq1name+"-"+seq2name
					pairseqnameR = seq2name+"-"+seq1name
					if (pairseqname not in pairs) and (pairseqnameR not in pairs):
						outseqs = [seq1,seq2]
						outfile = "./indv_gene_alns/"+"_".join([gene,seq1name,seq2name])+".fasta"
						temp_out = open(outfile,"w")
						SeqIO.write(outseqs, outfile, "fasta")
						temp_out.close()
			   			pairs.append(pairseqname)
						pairs.append(pairseqnameR)

core_genes = get_ids(sys.argv[1])
make_alns(core_genes)

