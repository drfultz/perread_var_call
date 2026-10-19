#!/usr/bin/env python3
import sys
import csv


### Takes a paf alignment file (with AS tag, as from minimap2) and outputs the name of the highest scoring reference
### and the best AS alignment score tag while also filtering out low quality matches
### Built for use with 'minimap2 -c' which outputs paf with cigar and AS tags

paf_file = snakemake.input[0]
feature = snakemake.wildcards.feature
calls_output = snakemake.output[0]
variant_len_cov_min = snakemake.params.VARIANT_LEN_COV_MIN
variant_len_cov_best = snakemake.params.VARIANT_LEN_COV_BEST


def parse_score(AS_field):
	#This parses the AS tag field into just the value
	return int(AS_field.rsplit('AS:i:', 1)[1])
def calc_ref_covered(line):
	#calculate the (end_ref-mapping - start_ref-mapping)/ref_length
	return (float(line[8]) - float(line[7]))/float(line[6])
def classify_and_output(output, readname, variant, score, variant_len_covered, variant_len_cov_min, variant_len_cov_best):
	read_id, unit, start = readname.rsplit(":", 2)

	#determine output type and output
	if variant_len_covered > variant_len_cov_best:
		output_line = f"{read_id}\t{unit}\t{start}\t{feature}:{variant}\t{score}"
	elif variant_len_covered > variant_len_cov_min:
		output_line = f"{read_id}\t{unit}\t{start}\t{feature}:{variant}%\t{score}"
	else:
		output_line = f"{read_id}\t{unit}\t{start}\t{feature}:Unmapped\t{score}"

	output.write(output_line + "\n")


##### Keep the highest alignment score for each read and output the readname, reference (which variant), and score
with open(calls_output,"w") as output, open(paf_file,'r') as paf:
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)

	##### Collect the best scoring match for each readname
	best_hits = {}
	for row in paf_reader:
		readname = row[0]
		### If this is the first line for this readname, store its values
		if readname not in best_hits:
			if row[4] == "*":
				best_hits[readname] = {"variant": row[5], "score": "NA", "ref_length_covered": 0.0}
			else:
				best_hits[readname] = {"variant": row[5], "score": parse_score(row[14]), "ref_length_covered": calc_ref_covered(row)}
		### If the readname has been seen, compare the current line's score with the previous highest score.
		else:
			current_score = parse_score(row[14])
			if current_score > best_hits[readname]["score"]:
				best_hits[readname] = {"variant": row[5], "score": current_score, "ref_length_covered": calc_ref_covered(row)}

	##### Output the best scoring match for each readname
	for readname, hit in best_hits.items():
		classify_and_output(output, readname, hit["variant"], hit["score"], hit["ref_length_covered"], variant_len_cov_min, variant_len_cov_best)