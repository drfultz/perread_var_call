#!/usr/bin/env python3
import sys
import csv
import random

###

if len(sys.argv)!=4:
	print("\tERROR: please input a subject .feat file, a query .feat file, and a window size", file=sys.stderr)
	sys.exit(2)

feat_file1 = sys.argv[1]
feat_file2 = sys.argv[2]
SVkmer_win = int(sys.argv[3]) #Perhaps need to catch errors where window is larger than the input feat lines

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
with open(feat_file1, 'r') as feat_obj1, open(feat_file2, 'r') as feat_obj2:
	feature_lines1 = [[x[0],x[1].split(",")] for x in csv.reader(feat_obj1, delimiter = "\t")]
	feature_lines2 = [[x[0],x[1].split(",")] for x in csv.reader(feat_obj2, delimiter = "\t")]
	#print(feature_lines)
	###Itereate through each line


for line1 in feature_lines1:
	###Iterate through windows of line1
	for i in range(len(line1[1][:-SVkmer_win+1])):
		output_list = [line1[0]]

		###against each line including self
		for line2 in feature_lines2:
			#seperate each sequence
			for x in range(10):
				output_list.append("NA")

			###Iterate through windows of line2
			for j in range(len(line2[1][:-SVkmer_win+1])):
				#print(f"{line1[1][i:i+SVkmer_win]}\t{line2[1][j:j+SVkmer_win]}")
				if line1[1][i:i+SVkmer_win] == line2[1][j:j+SVkmer_win]:
					output_list.append("1")
				else:
					output_list.append("0")

		print("\t".join(output_list))

	for x in range(10):
		print("\t".join(["NA" for x in output_list]))
