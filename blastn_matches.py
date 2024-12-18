#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 12:54:07 2024

@author: alba.marino
"""

# Get good and bad matches from library comparison with reciprocal masking with get_family_summary_paper.sh
# Output directory "summary_files" of get_family_summary_paper.sh renamed e.g. as summary_files_aedes_autoref

# What the script does: take a tsv file with the family IDs of one category (fragment, present, good, perfect) and the corresponding library,
# look at the summary file of the consensus and check if the reference ID has the same TE order as the match ID
# if more than one match ID is present and TE orders are different, choose the one with the longer hit.

# example run: python3 blastn_matches.py aedes_autoref_NEW.tsv ~/Desktop/TEannotation_benchmarking/libraries/aedes_libs/aedes_autolib_edit_NEW.fa summary_files_aedes_autoref aedes_autoref_blastn_matches.tsv


import sys, csv, subprocess
import os.path
import pandas as pd
from io import StringIO
from Bio import SeqIO
import numpy as np

input_tsv = os.path.abspath(sys.argv[1])
input_lib = os.path.abspath(sys.argv[2])
inout_dir = os.path.dirname(input_tsv)
summaries_path=os.path.abspath(sys.argv[3]) 
output_file= os.path.abspath(sys.argv[4])


with open (input_tsv, "r"):
	df_input = pd.read_csv(input_tsv, sep='\t')
	df_input= df_input.reset_index() 


dict_matches = {}
class1_tosearch = '|'.join(["ClassI/", "LINE", "Penelope", "SINE", "LTR"])
class2_tosearch = '|'.join(["ClassII/", "DNA", "MITE", "RC"])

for index, row in df_input.iterrows():
	cons_id = row['consensus_id']
	ref_order = row['order']
	ovle = row['overlap_level']
	if ovle == "Absent":
		dict_matches[cons_id] = {'match': np.nan}
	else:
		with open (f"{summaries_path}/{cons_id}.summary", "r") as summfile:
			df_summ = pd.read_csv(summfile, sep='\t', names=['cons_id', 'cons_len', 'beg_ref', 'end_ref', 'query_id', 'div'])
	
			if (ref_order == "ClassI_Others" and df_summ['query_id'].str.contains(class1_tosearch).any()) or (ref_order in '|'.join(["ClassI_Others", "LINE", "Penelope", "SINE", "LTR"]) and df_summ['query_id'].str.contains("ClassI/").any()):
	
				dict_matches[cons_id] = {'match': 'Good'} # wobbly ClassI correspondence
	
			elif (ref_order == "ClassII_Others" and df_summ['query_id'].str.contains(class2_tosearch).any()) or (ref_order in '|'.join(["ClassII_Others", "DNA", "MITE", "RC"]) and df_summ['query_id'].str.contains("ClassII/").any()):
	
				dict_matches[cons_id] = {'match': 'Good'} # wobbly ClassII correspondence
				
			elif ref_order in '|'.join(["Penelope", "LINE"]) and df_summ['query_id'].str.contains('|'.join(["Penelope", "LINE"])).any():
	
				dict_matches[cons_id] = {'match': 'Good'} # wobbly LINE-Penelope correspondence
	
			elif df_summ['query_id'].str.contains(ref_order).any():
	
				dict_matches[cons_id] = {'match': 'Good'} # exact match is necessary in all other cases
	
			else:
	
				dict_matches[cons_id] = {'match': 'Bad'}

df_matches = pd.DataFrame.from_dict(dict_matches, orient='index')
df_matches = df_matches.reset_index().rename(columns={"index":"consensus_id"})
df_output = df_input.merge(df_matches, how='left').drop('index', axis=1)

df_output.to_csv(output_file, sep="\t", index=False)
