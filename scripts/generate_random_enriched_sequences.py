#!/usr/bin/env python3
import os, sys
import subprocess
from Bio import SeqIO
import numpy as np
import json
import argparse

parser = argparse.ArgumentParser(description="Generate random enriched DNA fragments using BLAST output", version="0.01")

parser.add_argument("-m", action="store", dest="blastmap", help="BLAST output result file")
parser.add_argument("-r", action="store", dest="reference", help="Reference seqeunces to get DNA fragment")
parser.add_argument("-f", action="store", dest="fragments", help="List of randome fragment lengths")
parser.add_argument("-n", action="store", dest="total", type=int,  default = 10,  help="Total number of fragments to generate")

options=parser.parse_args()

def read_numbers_from_file(filename):
	'''
	read list of number, one in each line and returns an array of numbers
	'''
	numbers = []
	with open(filename) as input:
		for line in input:
			line = line.rstrip()
			if line == "": continue
			numbers.append(int(line))
	return numbers

def get_selected_chromosomes(chromosomes):
	'''
	return selected chromosomes from reference
	'''
	selected_chromosomes = {}
	with open(options.reference, 'r') as ref:
		for record in SeqIO.parse(ref, 'fasta'):
			if record.id in chromosomes:
				selected_chromosomes[record.id]=str(record.seq)

	return selected_chromosomes

def get_bait_map_data(blastdata):
	'''
	return blast out data in dictionary

	'''
	data ={}
	with open(options.blastmap, 'r') as blast:
		for line in blast:
			line = line.rstrip()
			if line == "":
				continue
			else:
				linearray = line.split()
				reference = linearray[1]
				mapStart = linearray[8]
				mapEnd  = linearray[9]
				if reference in data.keys():
					data[reference]["mapStart"].append(int(mapStart))
					data[reference]["mapEnd"].append(int(mapEnd))
				else:
					data[reference] = {"mapStart":[int(mapStart)], "mapEnd":[int(mapEnd)]}
	return data.keys(), data


def get_forward_fragment(chromosomeseq, start, fragmentLength):
	'''	returns forward dna fragments'''

	return chromosomeseq[start: start + fragmentLength]

def get_reverse_fragment(chromosomeseq, start, fragmentLength):
	''' return reverse dna fragment'''
	return chromosomeseq[start - fragmentLength:start]

def extract_dna_fragments(chromosomeseq, start, strand, fragmentLength):
	'''
	extract subsequence from chromosomes in Reference
	'''
	if strand == 0: ## forward strand
		if len(chromosomeseq) - start > fragmentLength:
			return get_forward_fragment(chromosomeseq, start, fragmentLength)
		else:
			return get_reverse_fragment(chromosomeseq, start, fragmentLength)

	elif strand == 1:
		if start > fragmentLength :
			return get_reverse_fragment(chromosomeseq, start, fragmentLength)
		else:
			return get_forward_fragment(chromosomeseq, start, fragmentLength)


def main():
	#get blast mapping data
	chromosomes, mapdata = get_bait_map_data(options.blastmap)
	# get the chromosomes with enriched DNA fragments
	selected_chromosomes = get_selected_chromosomes(chromosomes)

	random_fragment_lengths = read_numbers_from_file(options.fragments)

	fragment_counter = 0
	for fragmentLength in random_fragment_lengths:

		randomFragmentCreated=False
		while not randomFragmentCreated:
			random = np.random.random_integers(0, len(chromosomes)-1)
			chr_selected = chromosomes[random]

			number_of_mapped_regions_in_selected_chromosome = len(mapdata[chr_selected]['mapStart'])

			# choose a mapped region randomly
			if number_of_mapped_regions_in_selected_chromosome > 1:
				random_region = np.random.random_integers(0, number_of_mapped_regions_in_selected_chromosome-1)
			else:
				random_region = 0	# if there is only one region, it is in th 0th position in the array


			# choose randome start point in the mapped region
			start = mapdata[chr_selected]['mapStart'][random_region]
			end = mapdata[chr_selected]['mapEnd'][random_region]
			if start > end :
				start, end = end, start
			#print (start, end)
			#print("chromosome length ", len(selected_chromosomes[chr_selected]), " to get fragment size ", fragmentLength)
			random_start_position = np.random.random_integers(start, end)

			#choose forward or reverse strand, randomely, 0 - forward, 1 -reverse
			strand = np.random.random_integers(0,1)

			# now extract subsequence from reference
			fragment = extract_dna_fragments(selected_chromosomes[chr_selected], start, strand, fragmentLength)
			if len(fragment) < fragmentLength:
				continue
			else:
				print('>fragment' + str(fragment_counter) + " Length=" + str(fragmentLength))
				print(fragment)
				randomFragmentCreated = True
				fragment_counter += 1


if __name__ == '__main__':
	main()

exit(0)
