#!/usr/bin/env python3
import sys
import re

#aggregate multiple var call outputs into a features format with unit start positions

#merged_calls_file = sys.argv[1]
#features_file = "%s.feat" % merged_calls_file.rsplit('.', 1)[0]

merged_calls_file = snakemake.input[0]
features_file = snakemake.output[0]

with open(merged_calls_file,'r') as merged_calls, open(features_file,"w") as features_out:
	working_readname = f"{merged_calls_file}" # makes the first output line start with the filename
	working_features_string = "feature_variant_calls" # column name in first line
	working_startpos_string = "unit_start_pos" # column name in first line
	for input_line in merged_calls:
		#print(repr(line))
		#current_readname_match = re.match(r"[^\t]+", line) #every character until first tab
		current_readname = input_line.split("\t")[0]
		#print(current_readname)
		if current_readname == working_readname:
			next_features = re.findall(r"[^\t]+:[^\t\n]+", input_line)
			next_string = "&".join(next_features)
			#print(next_string)
			working_features_string = working_features_string + "," + next_string
			#print(working_features_string)
			working_startpos_string = working_startpos_string + "," + input_line.split("\t")[2]
		else:
			#write the finished line
			output_line = f"{working_readname}\t{working_features_string}\t{working_startpos_string}"
			print(output_line)
			features_out.write(output_line + "\n")
			#reset
			working_readname = current_readname
			#append
			features_list = re.findall(r"[^\t]+:[^\t\n]+", input_line) #look for but exclude tabs on first side, tab or newline on second side
			working_features_string = "&".join(features_list)
			#print(current_readname)
			#print(working_features_string)
			working_startpos_string = input_line.split("\t")[2]
	print(f"{working_readname}\t{working_features_string}\t{working_startpos_string}")
	features_out.write(f"{working_readname}\t{working_features_string}\t{working_startpos_string}" + "\n")