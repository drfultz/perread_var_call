#!/usr/bin/env python3
import sys
import csv

###

if len(sys.argv)!=3:
	print("\tERROR: please input a .feat file and a window size", file=sys.stderr)
	sys.exit(2)

feat_file = sys.argv[1]
SVkmer_window = int(sys.argv[2])

###Check if window size is odd
if SVkmer_window % 2 == 0:
	print("\tERROR: window size must be odd", file=sys.stderr)
	sys.exit(2)
else:
	win_reach = (SVkmer_window-1)//2

###read in the feat file
with open(feat_file, 'r') as feat_obj:
	for line in csv.reader(feat_obj, delimiter = "\t"):
		print(line[0])

		for x in range(win_reach):
			print("NA")

		units = line[1].split(",")
		for i,unit in enumerate(units[win_reach:(win_reach*-1)], start=win_reach):
			#print(unit.split("&"))

			hom_score = 1
			for r in range(1,win_reach+1):
				if unit == units[i-r]:
					hom_score += 1
				if unit == units[i+r]:
					hom_score +=1
			print(hom_score)

		for x in range(win_reach):
			print("NA")
