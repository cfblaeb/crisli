"""
ok heres the plan
This is a tool for a creating gRNA libraries. First I will focus on crispr inhibition libraries
then later we I may extend it if needed and any generalities I discover will be dealt with then, not now.
That means I wont even consider the other potential future uses right now.
It means that if you want to do anything other than generate a crispr i library for bacteria, then carefully go over the source
Especially the filtering stuff should probably be generalized by putting it into functions and then maybe via options the user could select a filter/selection package.

This is going to be a command line tool using standard python command line tool library (argparse)

As input it will need (for inhibition library)
external inputs:
1) A sequence file (fasta)
2) A list of regions in which to find gRNAs
	The crisli region file format is a tab separated file with 1 region per line (I suppose this could have been a standard annotation format...)
	id<tab>contig name<tab>strand<tab>start position<tab>end position
	The strand indicates where to look for gRNA
	1 = search on forward strand
	0 = search both strands
	-1 = search reverse strand
3) gRNAs to find per region

internal inputs: (input that maybe should be external options)
1) Filters: GC (20-80%), selection strategy


# search strategy (this is based on the regions being Transcription Units (TU) and greatest effect is seen early in the region
# and we dont care much that a grna targets sequence outside of the defined regions
	find all grnas in regions and associate grnas with their region (a grna may reside in overlapping regions)
	remove identical/similar ones (same sequence but different position) so only unique ones are left (grna may still be identical to sequence outside the regions)
	score remaining with crispy++
	pick X best grnas per region based on position+score or if less than X available, pick all
	specifically:
		filter based on 20% < GC < 80%
		filter based on crispy score > 99th percentile
		if still more than args.n left, then pick the earliest in the region

so, there we are, a somewhat defined idea
"""

import argparse
import sqlite3
from subprocess import run
from shlex import split as shsplit
from Bio import SeqIO
from helpers.filters import gc_filter, offtarget_score_filter


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
	parser.add_argument("n", type=int, default=4, help="Number of gRNA per region")
	parser.add_argument("s", type=int, default=0, help="search strategy. 0: if strand is -1 prefer near region start, if strand is 1 prefer near region end")
	parser.add_argument("seq_file", type=argparse.FileType(), help="fasta file")
	parser.add_argument("region_file", type=argparse.FileType(), help="file with regions to search for gRNA")

	args = parser.parse_args()

	# additional settings

	# size of gRNA
	sog = 20  # this effects crispy++ scoring...theres a hardcoded size in there

	print("read commandline. Reading fasta.")
	seqrecs = SeqIO.to_dict(SeqIO.parse(args.seq_file, "fasta"))
	print("Read fasta.")
	args.seq_file.close()
	print("Creating results.db (sqlite database, will overwrite existing)")
	dbc = sqlite3.connect('results.db')
	# create table structure
	dbc.executescript('''
		DROP TABLE IF EXISTS contigs;
		DROP TABLE IF EXISTS region;
		DROP TABLE IF EXISTS grna;
		DROP TABLE IF EXISTS offscores;
		CREATE TABLE contigs (
			id TEXT PRIMARY KEY,
			seq TEXT);
		CREATE TABLE region (
			id TEXT PRIMARY KEY,
			contig TEXT,
			strand TEXT,
			spos int,
			epos int,
			FOREIGN KEY(contig) REFERENCES contigs(id));
		CREATE TABLE grna (
			id INTEGER PRIMARY KEY,
			seq text,
			strand text,
			pos int,
			region TEXT,
			FOREIGN KEY(region) REFERENCES region(id));
		CREATE TABLE offscores (
			id INTEGER PRIMARY KEY,
			offscore int,
			grna INTEGER,
			FOREIGN KEY (grna) references grna(id));
	''')

	print("looking for grnas in regions from region file")
	# go over tss's and find grnas around them that fulfill some requirements
	contigs = set()  # store contigs
	regions = set()	 # store tss's
	grnas = set()	 # store grnas
	with args.region_file as f:
		for line in f.readlines():  # each line is a tss
			region_id, contig, strand, spos, epos = line.strip().split("\t")
			spos = int(spos)
			epos = int(epos)
			contigs.add((contig, str(seqrecs[contig].seq), ))
			regions.add((region_id, contig, strand, spos, epos))

			# find each grna in the region
			region_seq = seqrecs[contig][spos:epos].seq
			if strand == '1':  # search the forward strand
				for pot_pos in findall("GG", region_seq):
					# get grna seq without pam
					grna_5mark_pos = spos+pot_pos-sog-1
					gseq = seqrecs[contig].seq[grna_5mark_pos:grna_5mark_pos+sog]
					grnas.add((str(gseq).upper(), "1", grna_5mark_pos, region_id))
			else:  # reverse strand
				# find CC's and store them reverse transcribed
				for pot_pos in findall("CC", region_seq):
					# get grna seq without pam
					grna_5mark_pos = spos + pot_pos + 3 + sog
					gseq = seqrecs[contig].seq[grna_5mark_pos-sog:grna_5mark_pos].reverse_complement()
					grnas.add((str(gseq).upper(), "-1", grna_5mark_pos, region_id))

	print("writing results to db")
	# write contig, tss and grnas to db
	dbc.executemany("INSERT INTO contigs VALUES (?, ?)", contigs)
	dbc.executemany("INSERT INTO region VALUES (?, ?, ?, ?, ?)", regions)
	dbc.executemany("INSERT INTO grna (seq, strand, pos, region) VALUES (?, ?, ?, ?)", grnas)
	dbc.commit()

	print("Removing duplicate gRNA that are present at multiple locations (gRNA may also exists as duplicates due to overlapping regions, these will be allowed to persist)")
	dbc.execute("delete from grna where seq in (select s.seq from grna s inner join grna t on s.seq = t.seq and s.pos != t.pos group by s.seq)")

	print("Write remaining gRNA to file for crispy++ scoring")
	with open('scoreme', 'w') as f:
		f.writelines(f"{grna}\t{gid}\n" for gid, grna in dbc.execute("select id, seq from grna"))

	print("Score with crispy++. Note that if you are running this on a network drive its gonna run slow.")
	run(shsplit("crispy score 8 scoreme"))
	print("reading scores")
	with open('scored.laeb') as f:
		scored = [x.strip().split("\t") for x in f.readlines()]  # seq, gid, score

	print("storing in results.db")
	dbc.executemany("INSERT INTO offscores (offscore, grna) VALUES (?, ?)", ((score, gid) for seq, gid, score in scored))
	dbc.commit()
	print(f"going over regions and in cases where there are more than {args.n} grnas a gc filter, the score and the position will be used to select")
	#first lets get all grnas....I think its ok to keep them in memory and faster than retrieving them individually
	region_grna = {}
	for gid, seq, strand, pos, region, offscore in dbc.execute(("select g.id, g.seq, g.strand, g.pos, g.region, o.offscore from grna g join offscores o on g.id = o.grna")):
		if region in region_grna:
			region_grna[region].append([gid, seq, strand, pos, offscore])
		else:
			region_grna[region] = [[gid, seq, strand, pos, offscore]]

	# find regions with more than n grnas and decide which to keep
	# a threshold to use for offscore filtering
	threshold = dbc.execute("select offscore from offscores order by offscore asc limit 1 offset(select count(*) from offscores) * 99/100-1").fetchone()[0]  # 99th percentile
	for region_id, region_strand, grna_count in dbc.execute("select region.id, region.strand, COUNT(grna.id) c from region INNER JOIN grna on grna.region = region.id GROUP BY region.id having c > 4"):
		# try to remove grna using filters and if at any point theres only args.n left, then stop.
		# If after all filters there are still > args.n then pick the 5' most ones
		# filter 1: position must be in first half region
		# filter 2: %GC
		# filter 3: offscore
		# select args.n 5' most grnas left (strand dependent)

		grnas = region_grna[region_id]  # ok lets get the grnas

		filter_functions = [(gc_filter, {}), (offtarget_score_filter, {'threshold': threshold})]
		filtered_grna_ids = []
		for func, opts in filter_functions:
			grnas_left = [x for x in grnas if x[0] not in filtered_grna_ids]
			if len(grnas_left) > args.n:
				filtered_grna_ids += func(grnas, args.n, **opts)

		grnas_left = [x for x in grnas if x[0] not in filtered_grna_ids]

		if len(grnas_left) > args.n:  # if the remaining amount of grnas is still larger than n
			# then pick the args.n 5' most targets
			# that would be, if strand is 1 then its the highest pos and if strand is -1 then its the lowest pos
			# in the current version of grna from 1 region will have same strand
			re = False if grnas_left[0][2] == '1' else True  # determine region strand
			grnas_left = sorted(grnas_left, key=lambda x: x[3], reverse=re)  # sort by pos
			for gid, seq, strand, pos, offscore in grnas_left:
				filtered_grna_ids.append(gid)
				if len(grnas) - len(filtered_grna_ids) == args.n:  # if the remaining amount of grnas is now n
					break

		# ok end of filtering. There are now args.n grna
		if filtered_grna_ids:
			dbc.executemany("delete from grna where id is ?", [(x, ) for x in filtered_grna_ids])
	dbc.commit()

	print(f"Filtering complete. All regions now have {args.n} or less gRNAs.")
	print("preparing results from db for export")
	grnas = list(dbc.execute("select * from grna"))
	grna_dict = {}
	grna_to_region_dict = {}
	region_to_grna_dict = {}
	for id, seq, strand, pos, region in grnas:
		grna_dict[seq] = f'{strand}\t{pos}'
		if seq in grna_to_region_dict:
			grna_to_region_dict[seq] += f", {region}"
		else:
			grna_to_region_dict[seq] = region

		if region in region_to_grna_dict:
			region_to_grna_dict[region] += f", {seq}"
		else:
			region_to_grna_dict[region] = seq
	print("writing grna.tsv")
	with open('grna.tsv', 'w') as f:
		f.write("seq\tstrand\tpos\n")
		for key, val in grna_dict.items():
			f.write(f"{key}\t{val}\n")
	print("writing grna_to_region.tsv")
	with open('grna_to_region.tsv', 'w') as f:
		f.write("seq\tregion(s)\n")
		for key, val in grna_to_region_dict.items():
			f.write(f"{key}\t{val}\n")
	print("writing region_to_grna.tsv")
	with open('region_to_grna.tsv', 'w') as f:
		f.write("region\tgrna(s)\n")
		for key, val in region_to_grna_dict.items():
			f.write(f"{key}\t{val}\n")
	dbc.close()
