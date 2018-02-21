"""
ok heres the plan
This is a tool for a creating gRNA libraries. First I will focus on crispr inhibition libraries
then later we I may extend it if needed and any generalities I discover will be dealt with then, not now.
That means I wont even consider the other potential future uses right now.

This is going to be a command line tool using standard python command line tool library (argparse?)

As input it will need (for inhibition library)
external inputs:
1) A DNA sequence
2) A set of TSS with coordinates+strand for the DNA sequence
3) gRNA per TSS
internal inputs: (input that may become options)
4) rules about where to look for targets (e.g. +/- 50 bp from TSS, both strands/only 1
5) Filters: GC (20-80%), no duplicates

maybe it will run crispy++ scoring aswell and pick best, if so
6) strategy to pick best

so, there we are, a somewhat defined idea
"""

import argparse
from Bio import SeqIO


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Welcome to crisli. This tool can be used to generate a gRNA library.")
	parser.add_argument("n", type=int, default=4, help="Number of gRNA per location")
	parser.add_argument("seq_file", type=argparse.FileType(), help="genbank/fasta file with genome")
	parser.add_argument("tss_file", type=argparse.FileType(), help="file with tss")

	args = parser.parse_args()
	seqrec = SeqIO.parse(args.seq_file, "fasta")

	print([x for x in seqrec])
	args.seq_file.close()