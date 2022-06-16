#!/usr/bin/env python3
import sys
import csv
import random

###

if len(sys.argv)!=4:
	print("\tERROR: please input a .feat file, a minimum window size, and a maximum window size", file=sys.stderr)
	sys.exit(2)

feat_file = sys.argv[1]
SVkmer_minwin = int(sys.argv[2])
SVkmer_maxwin = int(sys.argv[3]) #Perhaps need to catch errors where maxwin*2 is larger than the input feat lines

###Functions
def calc_dr(parsed_feats, win_size):

	###reset variable for each size
	tandem_dup_count = 0
	total_count = 0

	###Iterate through units, ending with last full window
	for i in range(len(parsed_feats[:(win_size*-2+1)])):

		#print(f"{parsed_feats[i:i+win_size]}\t{parsed_feats[i+win_size:i+win_size+win_size]}")
		if parsed_feats[i:i+win_size] == parsed_feats[i+win_size:i+win_size+win_size]:
			#print("match")
			tandem_dup_count += 1
		
		total_count += 1

	return tandem_dup_count, total_count

###read in the feat file
with open(feat_file, 'r') as feat_obj:
	###Itereate through each line
	for line in csv.reader(feat_obj, delimiter = "\t"):

		###Output sequence name as header
		print(line[0])

		###Parse units for iteration of each unit in a line
		units = line[1].split(",")

		#
		random_counts = []

		###Iterate through window sizes
		for size in range(SVkmer_minwin,SVkmer_maxwin+1):
			values = calc_dr(units, size)
			print(f"win_{size}\t{values[0]/values[1]}")
			#
			random_counts.append([])

		###randomize 1000 times and print averages
		for x in range(1000):
			random.shuffle(units)
			###Iterate through window sizes
			for i,size in enumerate(range(SVkmer_minwin,SVkmer_maxwin+1)):
				values = calc_dr(units, size)
				random_counts[i].append(values[0]/values[1])

				#print(f"{line[0]} random {x}\t{values[0]/values[1]}")

		for i,x in enumerate(random_counts, start=1):
			print(f"rand_win_{i}\t{sum(x)/len(x)}")

		print("\n")
