#!/usr/bin/env python3
import sys
import re
import csv

#aggregate multiple var call outputs into a features format

merged_calls_file = snakemake.input[0]
features_file = snakemake.output[0]

with open(merged_calls_file,'r') as merged_calls, open(features_file,"w") as features_out:
	#calls_reader = csv.reader(merged_calls, delimiter='\t', quoting=csv.QUOTE_NONE)
	working_readname = "%s" % merged_calls_file
	features_string = "feature variant calls"
	for line in merged_calls:
		#print(repr(line))
		current_readname_match = re.match(r"[^\t]+", line) #every character until first tab
		current_readname = current_readname_match.group()
		#print(current_readname)
		if current_readname == working_readname:
			next_features = re.findall(r"[^\t]+:[^\t\n]+", line)
			next_string = ",".join(next_features)
			#print(next_string)
			features_string = features_string + "&" + next_string
			#print(features_string)
		else:
			#write old
			print(working_readname + "\t" + features_string)
			features_out.write(working_readname + "\t" + features_string + "\n")
			#reset
			working_readname = current_readname
			#append
			features_list = re.findall(r"[^\t]+:[^\t\n]+", line) #match look for but exclude tabs on first side, tabor newline on second
			features_string = ",".join(features_list)
			#print(current_readname)
			#print(features_string)
	print(working_readname + "\t" + features_string)
	features_out.write(working_readname + "\t" + features_string + "\n")
