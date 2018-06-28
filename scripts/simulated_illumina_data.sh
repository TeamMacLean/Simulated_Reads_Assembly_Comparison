#!/bin/bash 

~/art_bin_MountRainier/art_illumina --seqSys HS25 --in data/PGSC_DM_v3_scaffolds.fasta --len 125 --paired --minQ 10 --maxQ 40 --fcov 1000 -o simulated_data/HiSeq2500_nblrr_enrichment_simulated_from_potato_assembly.fastq --id HiSeq2500 --sdev 20 --mflen 250
