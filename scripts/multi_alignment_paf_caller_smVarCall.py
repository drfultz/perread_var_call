#!/usr/bin/env python3
import sys
import subprocess
import csv


### Takes a paf alignment file (with AS tag, as from minimap2) and outputs the name of the highest scoring reference
### and the best AS alignment score tag while also filtering out low quality matches


### Built for use with 'minimap2 -c' which outputs paf with cigar and AS tags

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
		ref_length_cov = 0.0
		with open(paf_file,'r') as paf2:
			paf_reader2 = csv.reader(paf2, delimiter='\t', quoting=csv.QUOTE_NONE)
			### save the values for the best scoring match
			for row in paf_reader2:
				if readname == row[0] and row[4] != '*':
					current_score = int(row[14].rsplit('AS:i:', 1)[1]) #### this parses the AS tag field into just the value
					current_ref_cov = (float(row[8]) - float(row[7]))/float(row[6])
					#print("%s\t%s\t%s" % (row[5],current_score,current_ref_cov))
					if current_score > score:
						variant = row[5]
						score = current_score
						ref_length_cov = current_ref_cov

		### determine output type
		if ref_length_cov > 0.95:
			print("%s\t%s\t%s:%s\t%s" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature,variant,score))
			output.write("%s\t%s\t%s:%s\t%s\n" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature,variant,score))
		elif ref_length_cov == 0.0:
			print("%s\t%s\t%s:Unmapped\tNA" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature))
			output.write("%s\t%s\t%s:Unmapped\tNA\n" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature))
		else:
			print("%s\t%s\t%s:%s^\t%s" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature,variant,score))
			output.write("%s\t%s\t%s:%s^\t%s\n" % (readname.rsplit(':')[0],readname.rsplit(':')[1],feature,variant,score))
			

subprocess.run("sort %s -o %s" % (calls_output,calls_output), shell=True)
