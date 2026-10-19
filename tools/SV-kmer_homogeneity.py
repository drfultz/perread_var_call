#!/usr/bin/env python3
import sys
import csv

##### Calculate how similar each feature call is to the units around it within a given window
##### and output (print) the scores as a vertical 'track'
##### Homogeneity score is number of exact feature calls at the whole unit level within the window
##### Intended for use on a variant called contig, not large read sets.


if len(sys.argv)!=3:
	print("\n\tERROR: please input a window size and a .feat file.  Writes to file named after input in WD.\n", file=sys.stderr)
	sys.exit(2)

SVkmer_window = int(sys.argv[1])
feat_file = sys.argv[2]
outname = feat_file.rsplit("/")[-1] + f".{SVkmer_window}win_homog.tsv"

###Check if window size is odd
if SVkmer_window % 2 == 0:
	print("\tERROR: window size must be odd", file=sys.stderr)
	sys.exit(2)
else:
	win_reach = (SVkmer_window-1)//2

##### Read in the feat file and set output
with open(feat_file, 'r') as feat_obj, open(outname, 'w') as output:
	##### For each line in the .feat file...
	for line in csv.reader(feat_obj, delimiter = "\t"):
		### Parse units and write header for each input line
		units = line[1].split(",")
		output.write(f"unit\thomogeneity\t{line[0]}\n")

		### Error message and skip if short line
		if len(units) < SVkmer_window:
			print(f"{line[0]} line is shorter than the window size", file=sys.stderr)
			continue			

		##### loop the units of the line
		for i, unit in enumerate(units):
			### if unit in usable interior, calculate score
			if win_reach <= i < len(units) - win_reach:
				### step out in both positive and negative from center unit to match
				hom_score = 1
				for r in range(1,win_reach+1):
					if unit == units[i-r]:
						hom_score += 1
					if unit == units[i+r]:
						hom_score +=1

			### else unit is at edges, report NA
			else:
				hom_score = "NA"

			### Write score
			output.write(f"{i+1}\t{hom_score}\n")
