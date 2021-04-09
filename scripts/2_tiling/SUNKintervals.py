#!/usr/bin/env python
from Bio import SeqIO # for reading fastq files
import argparse # for writing user-friendly command-line interfaces
import sys 
from collections import Counter #for counting # of times a chromosome's kmer is found in the read
import pandas as pd
import glob
import itertools
import csv

parser = argparse.ArgumentParser(description="")
parser.add_argument("-t", "--table", help="SUNKsharing.py outfile",  type=argparse.FileType('r'), default="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/leftCen8Read/positions_of_SUNKs_within_leftcen8readtrim450kb_and_other_reads_20mer.txt")
parser.add_argument("-o", "--out", help="output table (queryId\tSUNK_pair_id\tquery_dist\ttarget_dist)", type=argparse.FileType('w'), default="")
parser.add_argument("-n", "--num", help="number of shared SUNKs (queryId:num_of_shared_SUNKs)", type=argparse.FileType('w'), default="")
args = parser.parse_args()

# this is a function to determine similarity in distance between query and target SUNKs
# currently, it is set to 1% difference in distance and distance >20 (removes neighboring SUNKSs)
def is_similar(dist_query, dist_target):
    diff = abs(dist_query-dist_target)
    return (1.0*diff/(dist_target) <=0.01) and (diff>20)

#read in df
df = pd.read_table(args.table)

unique_query_reads = list(df['q_name'].unique())
query_df = df[df['q_name']==unique_query_reads[0]]

sunkid_position = []
for idx, row in query_df.iterrows():
    sunkid_position.append((row['kmer_color'], row['q_kmer_pos']))

sunkid_pairs = itertools.product(sunkid_position, sunkid_position)

target_df = df[['t_name', 'kmer_color', 't_kmer_pos']].drop_duplicates()

#map query reads, their unique pairs of SUNKs, and their corresponding distances
query_sunk_pair_mapping = {}
for query_read in unique_query_reads:
    query_df = df[df['q_name']==query_read]
    sunkid_position = []
    for idx, row in query_df.iterrows():
        sunkid_position.append((row['kmer_color'], row['q_kmer_pos']))

    sunkid_pairs = itertools.product(sunkid_position, sunkid_position)
    sunk_pair_dist = {}
    for pairs in sunkid_pairs:
        sunkid_a = pairs[0][0]
        sunkid_b = pairs[1][0]
        sunk_a_pos = pairs[0][1]
        sunk_b_pos = pairs[1][1]
        if sunkid_a != sunkid_b:
            sunk_pair_list = sorted([sunkid_a, sunkid_b])
            sunk_pair_id = '-'.join(map(str,sunk_pair_list))
            #str(sunkid_a) + '-' + str(sunkid_b)
            if sunk_pair_id not in sunk_pair_dist:
                sunk_pair_dist[sunk_pair_id] = set()
            sunk_pair_dist[sunk_pair_id].add(abs(sunk_a_pos - sunk_b_pos))
    
    query_sunk_pair_mapping[query_read] = sunk_pair_dist

#query_sunk_pair_mapping

#map target reads, their unique pairs of SUNKs, and their corresponding distances
target_sunkid_position = []
for idx, row in target_df.iterrows():
    target_sunkid_position.append((row['kmer_color'], row['t_kmer_pos']))

target_sunkid_pairs = itertools.product(target_sunkid_position, target_sunkid_position)
target_sunk_pair_dist = {}
for pairs in target_sunkid_pairs:
    sunkid_a = pairs[0][0]
    sunkid_b = pairs[1][0]
    sunk_a_pos = pairs[0][1]
    sunk_b_pos = pairs[1][1]
    if sunkid_a != sunkid_b:
        sunk_pair_list = sorted([sunkid_a, sunkid_b])
        sunk_pair_id = '-'.join(map(str,sunk_pair_list))
        #str(sunkid_a) + '-' + str(sunkid_b)
        if sunk_pair_id not in target_sunk_pair_dist:
            target_sunk_pair_dist[sunk_pair_id] = set()
        target_sunk_pair_dist[sunk_pair_id].add(abs(sunk_a_pos - sunk_b_pos))

#target_sunk_pair_dist

#want to print this as args.out

matching_query_results = []
with args.out as f:
	for query in query_sunk_pair_mapping:
		q_t_keys_intersect = set(query_sunk_pair_mapping[query].keys()).intersection(target_sunk_pair_dist.keys())
		for sunk_pair_intersect in q_t_keys_intersect:
			if is_similar(list(query_sunk_pair_mapping[query][sunk_pair_intersect])[0], 
					list(target_sunk_pair_dist[sunk_pair_intersect])[0]):	
					args.out.write(("{}\t{}\t{}\t{}\n".format(query, sunk_pair_intersect, query_sunk_pair_mapping[query][sunk_pair_intersect], target_sunk_pair_dist[sunk_pair_intersect])))
					matching_query_results.append((query, 
						sunk_pair_intersect,
						query_sunk_pair_mapping[query][sunk_pair_intersect],
						target_sunk_pair_dist[sunk_pair_intersect]))


#matching_query_results
#len(matching_query_results)

#
query_unique_matched_pair_count = {}
for query_kmer_pair in matching_query_results:
    if query_kmer_pair[0] not in query_unique_matched_pair_count:
        query_unique_matched_pair_count[query_kmer_pair[0]]=set()
    query_unique_matched_pair_count[query_kmer_pair[0]].add(query_kmer_pair[1])

#want to print this as args.num

#query_unique_matched_pair_count_final = {}
with args.num as f:
	for query_kmer_pair in query_unique_matched_pair_count:
		if len(query_unique_matched_pair_count[query_kmer_pair]) > 1:
			#query_unique_matched_pair_count_final[query_kmer_pair] = len(query_unique_matched_pair_count[query_kmer_pair])
			args.num.write("{}:{}\n".format(query_kmer_pair, len(query_unique_matched_pair_count[query_kmer_pair])))
    
#query_unique_matched_pair_count_final
