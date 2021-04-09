This repository contains scripts for the targeted assembly of complex genomic regions using singly unique nucleotide k-mers (SUNKs), or strings of length k that occur once in the genome. The scripts provided here use a SUNK library (i.e. generated with Jellyfish) and a set of ONT reads as input to barcode ONT reads and filter them based on shared SUNKs and their pairwise distances.

All scripts are designed to run on the Eichler lab cluster and will need to be modified to run on different clusters.
