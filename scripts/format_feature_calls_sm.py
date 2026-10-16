#!/usr/bin/env python3
import sys
from collections import defaultdict

### For each readname (first field) and read unit, find feature calls from every input file,
### then join these units together to form a features field (2nd field). Unit start positions are also joined
### into a list as the third field. In this way, the feature calling for each read is output as a single line in
### a '.feat' file. Does not assume that all readname/unit pairings have a call for each feature.

#merged_calls_file = sys.argv[1]
#features_file = "%s.feat" % merged_calls_file.rsplit('.', 1)[0]

calls_files = snakemake.input
features_file = snakemake.output[0]

##### Iterate the input files, and collect the calls from the same readname / unit pairs into a nested dict
reads = defaultdict(dict)
for file_name in calls_files:
	with open(file_name, 'r') as cur_calls:
		### Iterate through each line from a file 
		for line in cur_calls:
			## Parse the line into fields
			read_id, unit, start, call, score = line.rstrip("\n").split("\t")
			## If the unit of a readname is not in the nested dict, add it and the start site
			if unit not in reads[read_id]:
				reads[read_id][unit] = (start, [])
			## Add the call
			reads[read_id][unit][1].append(call)


##### Unpack the nested dict into output formatting and write
with open(features_file,"w") as features_out:
	### loop each read
	for read_id, units in reads.items():
		unit_calls = []
		unit_starts = []
 		### Aggregate the units of a read in order
		for unit in sorted(units, key=int):
			unit_calls.append("&".join(units[unit][1]))
			unit_starts.append(units[unit][0])
 		### Format and output
		features_string = ",".join(unit_calls)
		startpos_string = ",".join(unit_starts)
		features_out.write(f"{read_id}\t{features_string}\t{startpos_string}\n")