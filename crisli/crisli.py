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
	The crisli tss format is a tab separated file with:
	1 tss per line
	id<tab>contig name<tab>strand<tab>position
3) gRNA per TSS
internal inputs: (input that may become options)
4) rules about where to look for targets (e.g. +/- 50 bp from TSS, both strands/only 1
5) Filters: GC (20-80%), no duplicates

maybe it will also run crispy++ scoring and pick best, if so
6) strategy to pick best

so, there we are, a somewhat defined idea
"""

import argparse
import sqlite3
from Bio import SeqIO


def findall(sub, string):
	index = 0
	while index < len(string):
		index = string.upper().find(sub, index)
		if index == -1:
			break
		yield index
		index += 1


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Welcome to crisli. This tool can be used to generate a gRNA library.")
	parser.add_argument("n", type=int, default=4, help="Number of gRNA per location")
	parser.add_argument("seq_file", type=argparse.FileType(), help="fasta file with genome")
	parser.add_argument("tss_file", type=argparse.FileType(), help="file with tss")

	args = parser.parse_args()

	# additional settings
	# gRNA scan zone upstream/downstream
	gszu = 35
	gszd = 50
	# size of gRNA
	sog = 20  # this effects crispy++ scoring...theres a hardcoded size in there

	seqrecs = SeqIO.to_dict(SeqIO.parse(args.seq_file, "fasta"))
	args.seq_file.close()

	# strategy suggestion 1:
	# find all grnas in promoter zone (tss +/- gszd/gszu)
	# remove duplicates


	conn = sqlite3.connect('results.db')
	c = conn.cursor()

	c.execute("DROP TABLE IF EXISTS contigs")
	c.execute("DROP TABLE IF EXISTS tss")
	c.execute("DROP TABLE IF EXISTS grna")

	c.execute('''CREATE TABLE contigs (
				id TEXT PRIMARY KEY, 
				seq TEXT
				)''')

	c.execute('''CREATE TABLE tss (
				id TEXT PRIMARY KEY,
				contig TEXT,
				strand TEXT, 
				pos int,
				FOREIGN KEY(contig) REFERENCES contigs(id)
				)''')

	c.execute('''CREATE TABLE grna (
				id INTEGER PRIMARY KEY,
				seq text, 
				strand text, 
				pos int,
				tss TEXT,
				FOREIGN KEY(tss) REFERENCES tss(id)
				)''')

	with args.tss_file as f:
		for line in f.readlines():
			id, contig, strand, pos = line.strip().split("\t")
			pos = int(pos)
			# check if contig already exists
			if c.execute("SELECT Count() FROM contigs WHERE id is ?", (contig,)).fetchone()[0] == 0:
				c.execute("INSERT INTO contigs VALUES (?, NULL)", (contig,))
				conn.commit()
			c.execute("INSERT INTO tss VALUES (?, ?, ?, ?)", (id, contig, strand, pos))
			conn.commit()

			if strand == '1':  # forward strand
				promoseq = seqrecs[contig][pos-gszu:pos+gszd].seq
				# find CC's and store them reverse transcribed
				for pot_pos in findall("CC", promoseq):
					# get grna seq without pam
					grna_start_pos = pot_pos+3+pos-gszu
					gseq = seqrecs[contig].seq[grna_start_pos:grna_start_pos+sog].reverse_complement()
					c.execute("INSERT INTO grna (seq, strand, pos, tss) VALUES (?, ?, ?, ?)", (str(gseq), "-1", grna_start_pos, id))
			else:  # reverse strand
				promoseq = seqrecs[contig][pos-gszd:pos+gszu].seq
			break

	conn.commit()
	conn.close()
