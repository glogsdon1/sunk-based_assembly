import os
import tempfile
import numpy as np
import pandas as pd
import json
import re
import glob
from pprint import pprint
from Bio import SeqIO

#
# setup the env for each exacution 
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
snake_dir = SNAKEMAKE_DIR + "/"
shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_python3.cfg; " % SNAKEMAKE_DIR)
python3 = snake_dir + "env_python3.cfg"
python2 = snake_dir + "env_python2.cfg"

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()



lines = open("reads.list").readlines()
reads = []
for line in lines:
	reads.append(line.strip())
print("Number of reads: {}".format(len(reads)))
IDS = range(0, len(reads))

print(len(IDS))

wildcard_constraints:
	    ID="\d+"

localrules: all, final


def ReadsFromNum(wildcards):
	#rtn = []
	#for key in wildcards:
	#	#print(key)
	#	rtn.append(reads[int(key)])
	rtn = reads[ int( wildcards) ]
	return(rtn)

rule all:
	input: "final"




rule read_pos:
	input:
		fofn="reads.list",
	output:
		pos = "pos/{ID}.pos",
	params:
		cluster=" -l mfree=1G -l h_rt=10:00:00  "
	run:
		readID = ReadsFromNum(wildcards["ID"])
		print(readID)
		shell("source ~/.bashrc; grep {readID} reads_sharing_20mers.pos > {output.pos}" )


rule SUNK_sharing:
	input: 
		pos = "pos/{ID}.pos",
	output: 
		txt = "pos/positions_of_{ID}.txt",
	params:
		cluster=" -l mfree=2G -l h_rt=128:00:00 -pe serial 1 "
	threads: 1
	shell:
		"""
		source ~/modules.sh
		/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/SUNKsharing.py \
			 -1 parent_20.pos \
			 -2 {input.pos} \
			 -l lengths_of_reads_sharing_SUNKs_20mer.txt \
			 --out {output.txt}
		"""

rule SUNK_intervals:
	input:
		txt = "pos/positions_of_{ID}.txt",
	output:
		intervals = "pos/SUNKintervals_{ID}.tbl",
		numbers = "pos/SUNKnums_{ID}.txt",
	params:
		cluster=" -l mfree=4G -l h_rt=128:00:00 -pe serial 1 "
	threads: 1
	shell:
		"""
		source ~/modules.sh
		/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/SUNKintervals.py \
			-t {input.txt} \
			-o {output.intervals} \
			-n {output.numbers}
		"""

rule cat_tables:
	input: 
		txt = expand("pos/positions_of_{ID}.txt", ID=IDS),
		intervals = expand("pos/SUNKintervals_{ID}.tbl", ID=IDS),
		numbers = expand("pos/SUNKnums_{ID}.txt", ID=IDS),
	output:
		txt_sum = "summary/positions.txt",
		intervals_sum = "summary/SUNKintervals.tbl",
		numbers_sum = "summary/SUNKnums.txt"
	params:
		cluster=" -l mfree=4G -l h_rt=128:00:00 -pe serial 1 "
	threads: 1
	shell:
		"""
		cat {input.txt} | grep -v "t_name" > {output.txt_sum}
		cat {input.intervals} > {output.intervals_sum}
		cat {input.numbers} > {output.numbers_sum}
		"""

rule final:
	input:
		#numbers=expand("pos/SUNKnums_{ID}.txt", ID=IDS),
		#intervals=expand("pos/SUNKintervals_{ID}.tbl", ID=IDS),
		numbers_sum = "summary/SUNKnums.txt"
	output:
		"final"
	shell:
		"""
		touch {output}
		"""


'''
#!/bin/bash
mkdir -p reads aln results

NUM=$1
RUN=$2
FOFN=/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/newanalysis/chm13.fofn
REF=/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/newanalysis/data/Lit_plus_RM2_plus_addLit.bed.fasta
# get the NUMth line from fofn and use that bax file
BAX=$(awk '{ if(NR==$NUM) {print $0} }' $FOFN)
echo $BAX

# output
FASTQ=reads/$NUM.trimmed.fastq
OUT=aln/$NUM_$RUN.bam
OUT_filt_length=aln/$NUM_$RUN_filt_length.bam
OUT_filt_length_chrom=aln/$NUM_$RUN_filt_length_chrom.bam
RES=results/$NUM_$RUN.txt
RES_filt_length=results/$NUM_$RUN_filt_length.txt
RES_filt_length_chrom=results/$NUM_$RUN_filt_length_chrom.txt


if [ ! -f $FASTQ ]; then
	echo "Making the fastq from bax"
	pls2fasta -trimByRegion $BAX -fastq $FASTQ
fi

'''


