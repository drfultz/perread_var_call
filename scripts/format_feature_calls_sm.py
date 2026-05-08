#!/usr/bin/env python3
import sys

### For each readname (first field), find feature calls from every input file, combine the calls for each unit together,
### then join these units together to form a features field (2nd field). Unit start positions are also joined
### into a list as the third field.  In this way, the feature calling for each read is output as a single line in
### a '.feat' file.
###
### Assumes lines with the same readname are all adjacent and ordered by unit.

#merged_calls_file = sys.argv[1]
#features_file = "%s.feat" % merged_calls_file.rsplit('.', 1)[0]

calls_files = snakemake.input
features_file = snakemake.output[0]

# format and write output
def write_line():
	output_line = f"{working_readname}\t{working_features_string}\t{working_startpos_string}"
	#print(output_line)
	features_out.write(output_line + "\n")
# extract and join the variant calls from the list made from the input files
def find_features(line):
	features = []
	for calls in line:
		features.append(calls[3])
	return "&".join(features)

### Iterate the input files, and join the lines with the same line number from each file into lists
merged_calls = []
for file_name in calls_files:
	with open(file_name, 'r') as calls:
		if not merged_calls: # if empty (first file)
			for line in calls:
				merged_calls.append([line.split("\t")])
		else: #subsequent files
			for i,line in enumerate(calls):
				merged_calls[i].append(line.split("\t"))

### extract feature calls from the formatted inputs and output them as one line per readname
with open(features_file,"w") as features_out:
	working_readname = None

	for i, joined_line in enumerate(merged_calls):
		current_readname = joined_line[0][0]

		if current_readname == working_readname:
			#Append additional features and positions			
			next_features_string = find_features(joined_line)
			working_features_string = working_features_string + "," + next_features_string
			working_startpos_string = working_startpos_string + "," + joined_line[0][2]

		else:
			#write the finished line from the previous readname if not first line
			if i != 0:
				write_line()

			#reset with initial values
			working_readname = current_readname
			working_features_string = find_features(joined_line)
			working_startpos_string = joined_line[0][2]

	###write the final line
	write_line()
