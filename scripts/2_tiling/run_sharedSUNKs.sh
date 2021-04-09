#!/bin/bash
set -e          # terminate script if command exits with a nonzero exit status
set -u          # abort script if a variables value is unset
set -o pipefail # return nonzero exit status if any program in a pipe returns nonzero exit status


######### User-defined parameters
parentRead='13bda6ad-8b32-4001-bf43-54a3f13d5e84'
parentDir='/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo'
allReads='/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo/mapping_sunks/summary/reads_nodup.fa'
sharedSUNKsDir='/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo/mapping_sunks/summary/shared_SUNKs.tbl'
readsPosition='/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo/mapping_sunks/summary/reads.positions'
kmers='/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo/mapping_sunks/summary/chr8_asat_ONT.kmers'
cutoff='100'
readLength='100'

#make the directory
echo "Making the directory"
cd ${parentDir}
mkdir ${parentRead}

#get the parent read fasta
echo "Getting the parent read fasta"
cd ${parentRead}
/net/eichler/vol26/home/glogsdon/software/seqfilter/seqfilter -i ${allReads} -l <(printf "${parentRead}\n") -o ./${parentRead}.fa
source /net/eichler/vol26/home/glogsdon/modules.sh
samtools faidx ${parentRead}.fa
echo "Parent read retrieved"

#retrieving shared SUNK file
mkdir SUNKs
source /net/eichler/vol26/home/glogsdon/modules.sh
echo "Identifying SUNKs within read"
/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/KmersPos.py -k ${kmers} -fa ${parentRead}.fa -o SUNKs/${parentRead}.fa -t SUNKs/${parentRead}_20.tbl --pos SUNKs/${parentRead}_20.pos

echo "Number of reads sharing SUNKs: \
`grep -c "${parentRead}" ${sharedSUNKsDir}`"

echo "Generating data for histogram"
grep "${parentRead}" ${sharedSUNKsDir} | sort -rnk3,3 > SUNKs/${parentRead}_sharedSUNKs_20mer.tbl
awk '{print $3}' SUNKs/${parentRead}_sharedSUNKs_20mer.tbl > SUNKs/${parentRead}_sharedSUNKs_20mer.count
echo "Histogram data generated. Data is located here:
`readlink -f SUNKs/${parentRead}_sharedSUNKs_20mer.count`"

#run repeatmasker
mkdir repeatmasker
module load perl/5.14.2
module load RepeatMasker/3.3.0
RepeatMasker \
-species human \
-e wublast \
-dir ${parentDir}/${parentRead}/repeatmasker \
-pa 24 \
${parentDir}/${parentRead}/${parentRead}.fa
awk '{$16=""; print $0}' repeatmasker/${parentRead}.fa.out > repeatmasker/${parentRead}.fa.15.out 


#make the SUNK sudirectory
echo "Making the SUNK subdirectory"
cd ${parentDir}/${parentRead}/SUNKs
rm -rf ${cutoff}_or_more
mkdir ${cutoff}_or_more

#retrieve the read IDs of those with SUNKs greater than the cut-off
echo "Retrieving the reads"
awk '$3 >= 100' ${parentRead}_sharedSUNKs_20mer.tbl | awk '{print $1}'  > ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.tmp.txt
awk '$3 >= 100' ${parentRead}_sharedSUNKs_20mer.tbl | awk '{print $2}'  >> ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.tmp.txt
cat ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.tmp.txt | sort | uniq > ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.txt
echo "k-mer cut-off: ${cutoff}"
echo "read length cut-off: ${readLength}"
echo "Number of reads sharing SUNKs: `wc -l ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.txt`"

#determine the lengths of the reads
grep "${parentRead}" ${parentDir}/${parentRead}/${parentRead}.fa.fai | awk 'BEGIN{OFS="\t"}{print $1,$2}' > ${cutoff}_or_more/lengths_of_reads_sharing_SUNKs_20mer.txt
grep -f ${cutoff}_or_more/reads_sharing_20mers_with_${parentRead}.txt ${allReads}.fai | awk 'BEGIN{OFS="\t"}{print $1,$2}'>> ${cutoff}_or_more/lengths_of_reads_sharing_SUNKs_20mer.txt

#take only reads above above a certain size
awk '$2 >= 100000'  ${cutoff}_or_more/lengths_of_reads_sharing_SUNKs_20mer.txt > ${cutoff}_or_more/lengths_of_${readLength}kb_reads_sharing_SUNKs_20mer.txt
echo "Number of reads above ${readLength} kbp: `wc -l ${cutoff}_or_more/lengths_of_${readLength}kb_reads_sharing_SUNKs_20mer.txt`"

#filter the pos file
echo "Filtering the pos file"
awk '{print $1}' ${cutoff}_or_more/lengths_of_${readLength}kb_reads_sharing_SUNKs_20mer.txt > ${cutoff}_or_more/${readLength}kb_reads_sharing_SUNKs_20mer.txt
grep -f ${cutoff}_or_more/${readLength}kb_reads_sharing_SUNKs_20mer.txt ${readsPosition} > ${cutoff}_or_more/${readLength}kb_reads_sharing_20mers_with_${parentRead}.pos

#copy files into the snakemake directory
echo "Preparing the snakemake directory"
mkdir ${cutoff}_or_more/snakemake_for_SUNK_intervals
cp /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/mapping_SUNKs/tiling_across_centromere/demo/scripts/* ${cutoff}_or_more/snakemake_for_SUNK_intervals/

#change the files in the snakemake
cd ${parentDir}/${parentRead}/SUNKs/${cutoff}_or_more/
cp ${readLength}kb_reads_sharing_SUNKs_20mer.txt snakemake_for_SUNK_intervals/reads.list
grep -v "${parentRead}" snakemake_for_SUNK_intervals/reads.list > temp && mv temp snakemake_for_SUNK_intervals/reads.list
cp ../${parentRead}_20.pos snakemake_for_SUNK_intervals/parent_20.pos
cp ${readLength}kb_reads_sharing_20mers_with_${parentRead}.pos snakemake_for_SUNK_intervals/reads_sharing_20mers.pos
cp lengths_of_${readLength}kb_reads_sharing_SUNKs_20mer.txt snakemake_for_SUNK_intervals/lengths_of_reads_sharing_SUNKs_20mer.txt

#run the snakemake
cd snakemake_for_SUNK_intervals/
./runsnake.sh
cd summary
echo "Files are ready!"
readlink -f positions.txt
readlink -f SUNKintervals.tbl
cat SUNKnums.txt


