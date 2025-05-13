#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def is_file(filename):
    if not os.path.isfile(filename):
        raise argparse.ArgumentTypeError(f"{filename} is not a file")
    return filename

def get_args():
    parser = argparse.ArgumentParser(description='Remove stop codons from fasta alignment.')
    parser.add_argument("aln", help="Path to input FASTA alignment file", type=is_file)
    return parser.parse_args()

args = get_args()
in_file = args.aln
out_file = os.path.splitext(in_file)[0] + "_noStopCodons.fasta"

stop_codons = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}

# Read all sequences into memory
records = list(SeqIO.parse(in_file, "fasta"))
alignment_length = len(records[0].seq)

# Find positions of stop codons (same across all records)
stop_positions = []

for i in range(0, alignment_length, 3):
    codons = [str(record.seq[i:i+3]).upper() for record in records]
    if any(codon in stop_codons for codon in codons):
        stop_positions.append((i, i+3))

# Remove duplicates and sort in reverse to delete from end
stop_positions = sorted(set(stop_positions), reverse=True)

# Delete stop codons
for start, end in stop_positions:
    for record in records:
        seq_list = list(str(record.seq))
        del seq_list[start:end]
        record.seq = Seq("".join(seq_list))

# Write cleaned records
SeqIO.write(records, out_file, "fasta")
print(f"Output written to {out_file}")
