#!/usr/bin/env python3
import sys
import re

### Find var calls for every feature in a unit and merge them into a features entry,
### then merge each features entry into a features field (2nd field). Unit start positions are
### also merged into a list as the third field.  In this way, each read is output as
### a single line in a '.feat' file.
###
### Assumes a sorted input.

#merged_calls_file = sys.argv[1]
#features_file = "%s.feat" % merged_calls_file.rsplit('.', 1)[0]

merged_calls_file = snakemake.input[0]
features_file = snakemake.output[0]


def write_line():
	output_line = f"{working_readname}\t{working_features_string}\t{working_startpos_string}"
	#print(output_line)
	features_out.write(output_line + "\n")
def find_features(line):
	#look for but exclude tabs on first side, tab or newline on second side
	features = re.findall(r"[^\t]+:[^\t\n]+", line)
	return "&".join(features)

with open(merged_calls_file,'r') as merged_calls, open(features_file,"w") as features_out:
	working_readname = f"{merged_calls_file}" # makes the first output line start with the filename
	working_features_string = "feature_variant_calls" # column name in first line
	working_startpos_string = "unit_start_pos" # column name in first line

	for input_line in merged_calls:
		current_readname = input_line.split("\t")[0]

		if current_readname == working_readname:
			#Append additional features and positions			
			next_features_string = find_features(input_line)
			working_features_string = working_features_string + "," + next_features_string
			working_startpos_string = working_startpos_string + "," + input_line.split("\t")[2]

		else:
			#write the finished line from the previous readname
			write_line()
			#reset with initial values
			working_readname = current_readname
			working_features_string = find_features(input_line)
			working_startpos_string = input_line.split("\t")[2]

	###write the final line
	write_line()
