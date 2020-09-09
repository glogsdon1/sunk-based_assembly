#!/usr/bin/env python
from Bio import SeqIO # for reading fastq files
import argparse # for writing user-friendly command-line interfaces
import sys 
from collections import Counter #for counting # of times a chromosome's kmer is found in the read
import pandas as pd

one_read = "/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/leftCen8Read/leftcen8read_trim450kb_20.pos"
other_reads = "/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/20mer/summary/subset_reads_50SUNKsormore.pos"
lens = "/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr8_CCS/kmersinONTreads/leftCen8Read/lengths_of_reads_sharing_SUNKs_20mer.txt"

parser = argparse.ArgumentParser(description="")
parser.add_argument("-1", "--pos1", help="position table of SUNKs in one read",  type=argparse.FileType('r'), default=one_read)
parser.add_argument("-2", "--pos2", help="position table of SUNKs in all other reads",  type=argparse.FileType('r'), default=other_reads)
parser.add_argument("-l", "--len", help="lengths of reads other reads",  type=argparse.FileType('r'), default=lens)
parser.add_argument("-o", "--out", help="output table (readId1\treadId2\tnumOfSharedSUNKs)", type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()

# read position table into a dictionary
namedict = {}
kmer_color = {}
posdict = {}
lendict = {}

with args.len as f:
	for line in f:
		name, length = line.strip().split()
		lendict[name] = int(length)

def readfile(myfile, first=False):
	if(first):
		count = 0

	with myfile as f:
		for line in f:
			name, kmer, pos = line.strip().split()
			
			# add kmer to set of associated with read name
			if(name not in namedict):
				namedict[name] = set()
			namedict[name].add(kmer)
			
			# associate a read name and kmer with a position 
			posdict[(name, kmer)] = pos
			
			# add a color for the kmer if it is the first read (file)
			if(first):
				kmer_color[kmer] = count
				count += 1

	names = list(namedict.keys())
	return(names)

oneread = readfile(args.pos1, first=True)
assert len(oneread) == 1, "First argument can only have one read!!!"


# add header
form = 8*"{}\t" + "{}\n"
args.out.write(form.format("t_name", "t_len",
							"q_name", "q_len",
							"shared", 
							"kmer", "kmer_color",
							"t_kmer_pos", "q_kmer_pos"))

name1 = oneread[0]
other_names = readfile(args.pos2)
other_names.remove(name1)
for name2 in other_names:
	# kmers that intersect between read1 and read2
	inter = namedict[name1].intersection(namedict[name2])
	num = len(inter)
	if(num > 0):
		for inter_kmer in inter:
			args.out.write( form.format(
				name1, lendict[name1],
				name2, lendict[name2],
				num, 
				inter_kmer, 
				kmer_color[inter_kmer],
				posdict[(name1, inter_kmer)],
				posdict[(name2, inter_kmer)]
				) )
		

args.out.close()


