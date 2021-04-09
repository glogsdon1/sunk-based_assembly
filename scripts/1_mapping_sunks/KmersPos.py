#!/usr/bin/env python
from Bio import SeqIO # for reading fastq files
import argparse # for writing user-friendly command-line interfaces
import sys
from collections import Counter #for counting # of times a chromosome's kmer is found in the read
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("-k", "--kmers", help="input unique kmers",  type=argparse.FileType('r'), default="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr6_CCS/25/Illumina_all_count_CCS_chr6_25_91to124.kmer")
parser.add_argument("-fa", "--fasta", help="input fasta file", type=argparse.FileType('r'), default="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/kmerAnalysis/chr6_CCS/kmersinONTreads/25mer/test.fasta")
parser.add_argument("-c", "--count", help="min num of kmers in read", type=int, default=1)
parser.add_argument("-o", "--out", help="output fasta (with kmercount in description)", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument("-p", "--pos", help="output table (readid\tkmer\tpos)", type=argparse.FileType('w'), default="pos.tbl")
args = parser.parse_args()

# read all kmers into a list
kmers = set()
with args.kmers as f:
	for line in f:
		kmers.add(line.strip())
klen = 0
for idx, kmer in enumerate(kmers):
	if(idx == 0):
		klen = len(kmer)
	assert klen == len(kmer), "ERROR: unequal kmer sizes! {} {}".format(klen, len(kmer) )

# read in a fasta file
fasta = SeqIO.parse(args.fasta, "fasta")

# search each line in fasta file for kmer
outrec = []
outadj = {}
for idx, read in enumerate(fasta):
	seq = read.seq
	counter = 0 
	for i in range(len(seq) - len(kmer) + 1):
		kmer = seq[i:(i+klen)]
		if(str(kmer) in kmers):	
			counter += 1
			if(read.id not in outadj):
				outadj[read.id] = set()
			outadj[read.id].add( (kmer,i) )
		elif(str(kmer.reverse_complement()) in kmers):
			counter += 1
			if(read.id not in outadj):
				outadj[read.id] = set()
			outadj[read.id].add( (kmer.reverse_complement(),i) )
			

	# we need to write the read or not 
	print("read#:{}\tkmercount:{}".format(idx, counter), file=sys.stderr)
	if counter >= args.count: 
		read.description = read.description + "\tNumber of kmers:{}".format(counter)
		outrec.append(read)


SeqIO.write(outrec, args.out , "fasta")

# some code to write out adj to a file 
for key in outadj:
	if(len(outadj[key]) >= args.count):

		# this makes args.pos
		lines = ""
		for kmer, pos in outadj[key]:
			lines += "{}\t{}\t{}\n".format(key, kmer, pos)
		args.pos.write(lines)

name = args.pos.name
args.pos.close()
df = pd.read_table(name,  sep="\t", names=["readid", "kmer", "pos"])
df.sort_values(by=["readid", "pos"], inplace=True)
df.to_csv(name, index=False, sep="\t", header=False)



