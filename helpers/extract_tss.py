"""
Tool to convert TableS4.xlsx from DOI: 10.1126/science.1206848 into a CSV file with a crisli compatible format
The crisli tss format is a tab separated file with:
1 tss per line
id<tab>contig name<tab>strand<tab>position

in particular case, the excel file has strand and position
the genome has only 1 contig NC_000964.3
"""

import pandas as pd
import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Tool to convert tableS4.")
	parser.add_argument("tableS4", type=argparse.FileType(), help="tableS4 xlsx")

	args = parser.parse_args()
	file_name = args.tableS4.name
	args.tableS4.close()
	df = pd.read_excel(file_name)
	with open(f"{file_name[:-5]}.crisliTSS.tsv", 'w') as f:
		for i, r in df.loc[:, ['id', 'strand', 'beginTU']].iterrows():
			if r.beginTU != -1:
				f.write(f"{r.id}\tNC_000964.3\t{r.strand}\t{r.beginTU}\n")
