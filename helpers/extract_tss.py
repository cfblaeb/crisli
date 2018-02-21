"""
Tool to convert TableS4.xlsx from DOI: 10.1126/science.1206848 into a CSV file with a crisli compatible format
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
	print(df.head())
