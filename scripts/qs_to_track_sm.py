#!/usr/bin/env python3
import csv

### Takes a paf alignment file (with AS tag, as from minimap2)

feat_file = snakemake.input[0]
paf_file = snakemake.input[1]
sel_var_file = snakemake.input[2]
track_output = snakemake.output[0]
target_feat = snakemake.wildcards.feature


##### Open paf and loop through to sum counts of hits by feature mapped to
with open(paf_file, "r") as paf:
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	prev_name = ""
	read_counts_dict = {}
	pass_reads = 0
	for pline in paf_reader:
		new_name = pline[0]
		#if duplicated readname, skip or if matching residues / reference length is low (inaccurate mapping), skip
		if new_name == prev_name:# or int(pline[9])/int(pline[6]) < 0.98:
			print(f"{new_name} skipped for being dup of {prev_name} or low mapping % {int(pline[9])/int(pline[6])}")
			continue

		pass_reads +=1
		#if reference name is not in dictionary, create entry, else add one to entry
		if pline[5] not in read_counts_dict:
			read_counts_dict[pline[5]] = 1
		else:
			read_counts_dict[pline[5]] += 1

		#reset for next loop	
		prev_name = new_name

##### Open feat file for initial loop in order to sum counts of feat calls
with open(feat_file, "r") as feat:

	### first, find which field the target feature is in
	#read the first line and reset pointer
	first_line_list = feat.readline().split("\t")
	feat.seek(0)

	field_count = 0
	field_num = 0
	#loop through the call from first unit of the units field of the first line of the feat file
	for call in first_line_list[1].split(",")[0].split("&"):
		if call.split(":")[0] == target_feat:
			field_num = field_count
			break
		field_count += 1
	# no break
	else:
		print(f"\n\tERROR: the target feature {target_feat} was not found in the first set of calls in the feature_calls file.\n", file=sys.stderr)
		sys.exit(2)


	### build dict of counts of calls used in feature_calls file
	# open reader for feat file
	feat_reader1 = csv.reader(feat, delimiter='\t', quoting=csv.QUOTE_NONE)
	featcall_dict = {}
	pass_calls = 0
	# loop feature calls of every line
	for fline1 in feat_reader1:
		for calls1 in fline1[1].split(","):
			variant = calls1.split("&")[field_num].split(":")[1]
			#if value is Unmapped or ^, don't count in total or add to dict
			if variant == "Unmapped" or variant[-1] == "^":
				print(f"{variant} not used")
				continue

			pass_calls +=1

			#if call name is not in dictionary, create entry, else add one to entry
			if variant not in featcall_dict:
				featcall_dict[variant] = 1
			else:
				featcall_dict[variant] += 1


##### Calculate normalized values
norm_dict = {}
#for each feature call that was present in Chr/reference
for variantd in featcall_dict:
	#if there were reads mapped to this feature call
	if variantd in read_counts_dict:
		ratio = (read_counts_dict[variantd]/pass_reads)/(featcall_dict[variantd]/pass_calls)
		scaled_ratio = ratio/(ratio + 1)
		norm_dict[variantd] = scaled_ratio
	#else, nothing mapped gets a 0
	else:
		norm_dict[variantd] = 0


##### Process feat file and write outputs as it goes
with open(feat_file, "r") as feat, open(sel_var_file, "r") as sel_feats, open(track_output,"w") as output:
	# open reader for feat file
	feat_reader2 = csv.reader(feat, delimiter='\t', quoting=csv.QUOTE_NONE)
	sel_var_list = [x.rstrip() for x in sel_feats]
	# loop lines of feat file 
	for fline in feat_reader2:
		#each line (chromosome) gets it's own count reset (however, output merges all lines into unbroken columns) 
		count = 0
		#loop through the feature call fields of each line and output a normalized value for each match
		for calls,pos in zip(fline[1].split(","),fline[2].split(",")):
			count += 1

			name = f"{fline[0]}:{str(count).zfill(4)}"

			call_val = calls.split("&")[field_num].split(":")[1]
			if call_val in norm_dict:
				# positive values are part of select variants
				if call_val in sel_var_list:
					value = norm_dict[call_val]
				# calls not in select variants are given negative values so that data is not lost, but is easily ignored in R
				else:
					value = -(norm_dict[call_val])
			else:
				value = "NA"
			
			output.write(f"{name}\t{pos}\t{value}\n")
