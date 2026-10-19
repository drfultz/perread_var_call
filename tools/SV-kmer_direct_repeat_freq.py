#!/usr/bin/env python3
import sys
import csv
import random

### Takes a .feat file and, for each read, measures how often a block of consecutive units is immediately
### followed by an identical block (a direct tandem repeat of units).

### Inputs are expected to be 1) a .feat file from the feature-calling pipeline, 2) a minimum window size, and
### 3) a maximum window size, where window size is the number of units in each block being compared.
### The intended input is a feature-called contig, but it may work with reads.  Writes to stdout.

### For each read and window size w, the script slides along the read and compares each block of w units
### to the w units immediately after it. For each read, the frequency of observed direct repeats at each size
### is returned alongside the frequency when that read is randomly shuffled 1000X.

if len(sys.argv)!=4:
	print("\tERROR: please input a .feat file, a minimum window size, and a maximum window size", file=sys.stderr)
	sys.exit(2)

feat_file = sys.argv[1]
SVkmer_minwin = int(sys.argv[2])
SVkmer_maxwin = int(sys.argv[3]) #Perhaps need to catch errors where maxwin*2 is larger than the input feat lines


def is_simple_repeat(block):
	"""True if the block is not a shorter pattern repeated"""
	w = len(block)
	for p in range(1, w):
		if w % p != 0:
			continue
		pattern = block[:p]
		if pattern * (w // p) == block:
			return True
	return False

def calc_dr(parsed_feats, win_size):
	##### Calculate the number of direct tandem repeats of units

	### Set variable for each size
	tandem_dup_count = 0
	total_count = 0

	###Iterate through units, ending with last full window
	for i in range(len(parsed_feats[:(win_size*-2+1)])):
		###Check if the current window is made of simpler repeats and skip if so
		cur_win = parsed_feats[i:i+win_size]
		if is_simple_repeat(cur_win):
			continue
		### Look for direct tandem repeats
		if parsed_feats[i:i+win_size] == parsed_feats[i+win_size:i+win_size+win_size]:
			tandem_dup_count += 1
		
		total_count += 1
	### Return both the observed counts and the total possible windows
	return tandem_dup_count, total_count


###read in the feat file
with open(feat_file, 'r') as feat_obj:
	### Write a header
	print("read\twindow\tobserved\trandom_mean")
	sizes = range(SVkmer_minwin,SVkmer_maxwin+1)

	###Itereate through each line
	for line in csv.reader(feat_obj, delimiter = "\t"):

		###Parse units for iteration of each unit in a line
		units = line[1].split(",")

		observed = {}
		random_counts = {}

		###Iterate through window sizes and store observed values
		for size in sizes:
			dr_count, total = calc_dr(units, size)
			observed[size] = dr_count/total
			#Make empty list per size for later
			random_counts[size] = []

		###randomize 1000 times and store averages for each size
		for x in range(1000):
			random.shuffle(units)
			###Iterate through window sizes
			for size in sizes:
				dr_count, total = calc_dr(units, size)
				random_counts[size].append(dr_count/total)

		for size in sizes:
			mean_random = sum(random_counts[size])/len(random_counts[size])
			print(f"{line[0]}\t{size}\t{observed[size]}\t{mean_random}")
