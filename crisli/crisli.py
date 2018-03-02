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
	parser.add_argument("seq_file", type=argparse.FileType(), help="genbank file with genome sequence and annotation")
	parser.add_argument("tss_file", type=argparse.FileType(), help="file with tss")

	args = parser.parse_args()

	# additional settings
	# gRNA scan zone upstream/downstream
	gszu = 35
	gszd = 50
	# size of gRNA
	sog = 20  # this effects crispy++ scoring...theres a hardcoded size in there

	seqrecs = SeqIO.to_dict(SeqIO.parse(args.seq_file, "genbank"))
	args.seq_file.close()

	# strategy suggestion 1:
	# find all grnas in promoter zone (tss +/- gszd/gszu)
	# remove duplicates

	# create db
	dbc = sqlite3.connect('results.db')
	# create table structure
	dbc.executescript('''
		DROP TABLE IF EXISTS contigs;
		DROP TABLE IF EXISTS tss;
		DROP TABLE IF EXISTS grna;
		CREATE TABLE contigs (
			id TEXT PRIMARY KEY,
			seq TEXT);
		CREATE TABLE tss (
			id TEXT PRIMARY KEY,
			contig TEXT,
			strand TEXT,
			pos int,
			FOREIGN KEY(contig) REFERENCES contigs(id));
		CREATE TABLE grna (
			id INTEGER PRIMARY KEY,
			seq text,
			strand text,
			pos int,
			tss TEXT,
			FOREIGN KEY(tss) REFERENCES tss(id));
	''')

	# go over tss's and find grnas around them that fulfill some requirements
	contigs = set()  # store contigs
	tsss = set()	 # store tss's
	grnas = set()	 # store grnas
	with args.tss_file as f:
		for line in f.readlines():  # each line is a tss
			id, contig, strand, pos = line.strip().split("\t")
			pos = int(pos)
			contigs.add((contig, str(seqrecs[contig].seq), ))
			tsss.add((id, contig, strand, pos))

			# find each grna around the tss
			if strand == '1':  # forward strand
				promoseq = seqrecs[contig][pos-gszu:pos+gszd].seq
				# find CC's and store them reverse transcribed
				for pot_pos in findall("CC", promoseq):
					# get grna seq without pam
					grna_start_pos = pos-gszu+pot_pos+3
					gseq = seqrecs[contig].seq[grna_start_pos:grna_start_pos+sog].reverse_complement()
					grnas.add((str(gseq), "-1", grna_start_pos, id))
			else:  # reverse strand
				promoseq = seqrecs[contig][pos-gszd:pos+gszu].seq
				# find GG's and store them
				for pot_pos in findall("GG", promoseq):
					# get grna seq without pam
					grna_end_pos = pos-gszd+pot_pos-1
					gseq = seqrecs[contig].seq[grna_end_pos-sog:grna_end_pos]
					grnas.add((str(gseq), "1", grna_end_pos-sog, id))

	# write contig, tss and grnas to db
	dbc.executemany("INSERT INTO contigs VALUES (?, ?)", contigs)
	dbc.executemany("INSERT INTO tss VALUES (?, ?, ?, ?)", tsss)
	dbc.executemany("INSERT INTO grna (seq, strand, pos, tss) VALUES (?, ?, ?, ?)", grnas)
	dbc.commit()

	# Strategy forward:
	#  Save gRNAs (I guess they are already saved in the "grnas" set)
	#  Then remove duplicates from the gRNA table
	#  Then find all tss's that has less than the requested number of grnas
	#    Then go through those tss's and go downstream from gszd and find the needed extra grnas
	#      Check that those gRNAS dont exist in the original grnas set
	#      Search region up to first CDS stop codon

	#remove duplicate grnas
	dbc.execute("DELETE FROM grna WHERE seq in (SELECT seq FROM grna GROUP BY seq HAVING COUNT(seq) != 1)")

	#get grna per tss count
	gpt = [x for x in dbc.execute("select tss.id, COUNT(grna.id) c from tss INNER JOIN grna on tss.id = grna.tss GROUP BY tss.id")]
	print(gpt[:3])

	print([x for x in dbc.execute("SELECT COUNT(*) FROM grna")])

	dbc.close()
