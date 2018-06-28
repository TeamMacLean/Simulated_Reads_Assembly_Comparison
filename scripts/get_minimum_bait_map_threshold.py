#!/usr/bin/python3

import os, sys
import json

blastfile=sys.argv[1]
blasthandle=open(blastfile)

baithits=dict()
for line in blasthandle:
	line=line.rstrip()
	blastarray=line.split()
	baitname=blastarray[0]
	bait_align_percentage=float(blastarray[2])
	align_length=int(blastarray[3])
	align_gaps=int(blastarray[5])
	if align_length == 120 and align_gaps == 0:
		#print(bait_align_percentage, align_length, align_gaps)
		if baitname in baithits.keys():
			if baithits[baitname] < bait_align_percentage:
				#print("previous ", baithits[baitname], " now ", bait_align_percentage)
				baithits[baitname] = bait_align_percentage
				#print("New value ", baithits[baitname])
				#exit(0)
		else:
			baithits[baitname] = bait_align_percentage

with open("results/max_individual_bait_hits.json", 'w') as baithit:
	json.dump(baithits, baithit, sort_keys=True, indent=4, separators=(',',';'))

for key in baithits.keys():
	print(key, baithits[key])

exit(0)
