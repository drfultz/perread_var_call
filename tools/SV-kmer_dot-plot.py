#!/usr/bin/env python3
import sys
import csv

##### Makes a similarity matrix of feature calls (rather than DNA sequence) with the input window size.
##### Similarity is scored on a scale from 0 to the window size based on 
##### how many whole repeat units exactly match for SV calls.  Prints output in a tab seperated format.

##### Can be used on multi-line .feat files
##### Can do self vs self if the same file is provided twice.

if len(sys.argv)!=4:
	print("\tERROR: please input a subject .feat file, a query .feat file, and a window size", file=sys.stderr)
	sys.exit(2)

feat_file1 = sys.argv[1]
feat_file2 = sys.argv[2]
SVkmer_win = int(sys.argv[3]) #Perhaps need to catch errors where window is larger than the input feat lines

###read in the feature files
### does not trim out the '^' mark for uncertain calls. That should be useful for accurate sequence.
with open(feat_file1, 'r') as feat_obj1, open(feat_file2, 'r') as feat_obj2:
	feature_lines1 = [[x[0],x[1].split(",")] for x in csv.reader(feat_obj1, delimiter = "\t")]
	feature_lines2 = [[x[0],x[1].split(",")] for x in csv.reader(feat_obj2, delimiter = "\t")]

### for each feature call 'read'
for line1 in feature_lines1:
	###Iterate through windows of SV calls and output one line per window
	for i in range(len(line1[1][:-SVkmer_win+1])):
		#start each output list (output line) with the name of the 'read'
		output_list = [line1[0]]
		#define the current window of SV calls for the query 
		unitI = line1[1][i:i+SVkmer_win]

		###compare against each feature call 'read'
		for line2 in feature_lines2:
			#before each new 'target add spacer "NA" values that will create gray in R
			for x in range(10):
				output_list.append("NA")

			###Iterate through windows of the subject and append a score.
			for j in range(len(line2[1][:-SVkmer_win+1])):
				#define the current window of SV calls for the subject
				unitJ = line2[1][j:j+SVkmer_win]

				#initialize score for this set of windows
				score = 0
				### compare unit by unit
				for x, y in zip(unitI,unitJ):
					if x == y:
						score += 1
				#add score to output line
				output_list.append(str(score))

		### print the output line for each window
		print("\t".join(output_list))

	for x in range(10):
		print("\t".join(["NA" for x in output_list]))
