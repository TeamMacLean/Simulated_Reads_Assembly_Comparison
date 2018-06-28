## Introduction

This is a project on RenSeq enrichment sequence assembly comparison with different tools.

## Objective

Simulation of RenSeq enriched sequence data and assemble with different assembly tools.

## Methods

1) Crop out only the DNA fragment band of interest i.e. renseq enriched from the raw complete gel image.
2) Analyse the cropped gel image for the intensity of DNA fragment after renseq enrichment. Get the mean of the intensities and cumulative proportions of pixels across each row in the image.
3) Get length in base pair for difference in each row of intensities (from plot)
4) Generate any number (N) of random  values from 0 to 1. for e.g. N can be 100,000
5) Generate sequence lengths from the random values using the upper and lower DNA band values in gel and the number of rows of intensities analysed.
6) Blast the renseq baits (120bp) to a real pacbio renseq enriched data to get the best worst percentage hit. The baits have to be hit within a read with no gap. Plot the hit percentage and with manual inspection, find the bait hit percentage.
7) Randomly generate DNA fragments from a reference sequence (like sonication in the lab).
8) Blast and filter out only the fragments that get hit with percentage equal or higher than the hit percentage from manual inspective above.
9) simulate pacbio reads and illumina reads
10) do assembly
