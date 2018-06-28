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

Blast baits sequences to 4 more pacbio ROI consensus data. The more data I have, the more accurate threshold I get.

