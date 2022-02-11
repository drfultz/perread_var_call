#!/usr/bin/env python3
import sys
import subprocess
import csv


### Takes a paf alignment file (with AS tag, as from minimap2) and outputs the name of the highest scoring reference
### and the best AS alignment score tag while also filtering out low quality matches


### Built for use with 'minimap2 -c' which outputs paf with cigar and AS tags

# if len(sys.argv)!=3:
# 	print("ERROR:  please provide a paf with AS tags and a name, no more, no less\n")
# 	sys.exit(2)

# paf_file = sys.argv[1]
# feature = sys.argv[2]
# calls_output = "%s_calls.tsv" % paf_file.rsplit('.', 1)[0]

paf_file = snakemake.input[0]
feature = snakemake.wildcards.feature
calls_output = snakemake.output[0]


def parse_score(AS_field):
	#This parses the AS tag field into just the value
	return int(AS_field.rsplit('AS:i:', 1)[1])
def calc_ref_covered(line):
	#calculate the (end_ref-mapping - start_ref-mapping)/ref_length
	return (float(line[8]) - float(line[7]))/float(line[6])
def classify_and_output():
	#determine output type and output
	if ref_length_covered > 0.95:
		output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:{variant}\t{score}"
		#print(output_line)
		output.write(output_line + "\n")
	elif ref_length_covered > 0.8:
		output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:{variant}^\t{score}"
		#print(output_line)
		output.write(output_line + "\n")
	else:
		output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:Unmapped\t{score}"
		#print(output_line)
		output.write(output_line + "\n")

##### Keep the highest alignment score for each read and output the readname, reference (which variant), and score
with open(calls_output,"w") as output:
	with open(paf_file,'r') as paf:
		paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)

		readname = "readname:unit#:start_pos"
		variant = "variant"
		score = "score"
		ref_length_covered = 0.0

		##### output the values for the best scoring match of the previous readname when a new readname is encountered
		for i, row in enumerate(paf_reader):
			### If changed to a new readname, output the old and reset the variables
			if row[0] != readname:
				#ouput the old as long as not first line
				if i != 0:
					classify_and_output()

				###reset the variables
				readname = row[0]
				variant = row[5]
				if row[4] == "*":
					score = "NA"
					ref_length_covered = 0.0
				else:
					score = parse_score(row[14])
					ref_length_covered = calc_ref_covered(row)
			### If not a new readname, compare the current line's score with the previous highest score.
			else:
				#This block may be largely unnecessary - minimap2 appears to output the highest AS first.
				#Unless other comparisons are made between mappings for QC, comparing scores is superfluous.
				current_score = parse_score(row[14])
				#print("%s\t%s\t%s" % (row[5],current_score,current_ref_cov))
				if current_score > score:
					variant = row[5]
					score = current_score
					ref_length_covered = calc_ref_covered(row)

		##### Output the final best score (as there is no new readname to trigger output).
		classify_and_output()
