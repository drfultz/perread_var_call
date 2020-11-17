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

### uses a set to gather each unique readname
with open(paf_file,'r') as paf:

	readnames = set()
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	for row in paf_reader:
		readnames.add(row[0])
		#print(row[0])
	print("reads gathered\n")

### Keep the hightest alignment score for each read and output the readname, reference name, and score
with open(calls_output,"w") as output:
	for readname in readnames:
		#print(readname)
		variant = "Error"
		score = 0
		ref_length_covered = 0.0
		with open(paf_file,'r') as paf2:
			paf_reader2 = csv.reader(paf2, delimiter='\t', quoting=csv.QUOTE_NONE)
			### save the values for the best scoring match
			for row in paf_reader2:
				if readname == row[0] and row[4] != '*':
					current_score = int(row[14].rsplit('AS:i:', 1)[1]) #### this parses the AS tag field into just the value
					#print("%s\t%s\t%s" % (row[5],current_score,current_ref_cov))
					if current_score > score:
						variant = row[5]
						score = current_score
						ref_length_covered = (float(row[8]) - float(row[7]))/float(row[6])

		### determine output type
		if ref_length_covered > 0.95:
			output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:{variant}\t{score}"
			#print(output_line)
			output.write(output_line + "\n")
		elif ref_length_covered == 0.0:
			output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:Unmapped\tNA"
			print(output_line)
			output.write(output_line + "\n")
		else:
			output_line = f"{readname.split(':')[0]}\t{readname.split(':')[1]}\t{readname.split(':')[2]}\t{feature}:{variant}^\t{score}"
			#print(output_line)
			output.write(output_line + "\n")
			

subprocess.run("sort %s -o %s" % (calls_output,calls_output), shell=True)
