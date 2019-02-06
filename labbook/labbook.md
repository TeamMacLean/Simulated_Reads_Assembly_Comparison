## since the project started in March 2018

1) I have tried to simulate Illumina reads using ART.
2) Got gel image of DNA fragments captured by baits designed for RenSeq sequencing.
3) Image analysis to get intensities, plotting of the intensities
4) Assessing the bp length of each point from the plot. Generating random number to get random sequence lengths

## 14 June

Get the minimum best bait hit to genome.

Conditions are :

1) full length bait alignment (this ensures the bait mapping within the fragment)
2) no gaps in alignment

## 19 June 2018

A pacbio raw consensus (ROI) data enriched using renseq baits obtained from Kamil. The data are sequenced in a RSII smrtcell. A DNA fragment is sequenced repeated in circular fashion. The subreads (part between the adapters) are collapsed into one consensus sequence. Baits mapping is done on the consensus sequence to get the threshold as they represent closest to the real DNA sequence due to their consensus nature.

## 20 June 2018

Blast baits sequences to 4 more PacBio ROI consensus data. The more data I have, the more accurate threshold I get.

## July 2018



## August 2018

Random fragment lengths generated

## 30 Sept 2018

Script ready to generate DNA fragments. DNA fragments generated.

## Oct 10

Generating illumina HiSeq2500 150 pb paired end reads with default parameters.

## Jan 20 2019

Writing a new pipeline to add more variablity in the simulated data and assembly.

## Jan 24 2019

Draft Pipeline completed. Running the pipeline to simulated a dataset for a set of parameters. It also does assembly of the data and map the baits to the assembly.

## Jan 31 2019

Jobs are runnning. Assembly and analysis of few datasets completed. Getting the report for the completed datasets.

