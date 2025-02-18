#!/usr/bin/env python3
import csv
import sys

### Takes a feat file and prints feat calls that are exclusive to one chromosome and within
### a max window of user provided distance (in bp)

#!! Currently assumes only a single feature is called in feat file which is not ideal


if len(sys.argv)!=3:
	print("\n\tERROR: please input a .feat file (with only a single feature) and an integer of max distance.\n", file=sys.stderr)
	sys.exit(2)

##### INPUTS
feat_file = sys.argv[1]

try:
	distance_cutoff = int(sys.argv[2])
except:
	print("\n\tERROR: The second positional argument should be an integer.\n", file=sys.stderr)
	sys.exit(2)

##### Open feat file and make a dictionary for each chromosome with each variant as a key paired to it's positions in that Chr
with open(feat_file, "r") as feat:
	# open reader for feat file
	feat_reader1 = csv.reader(feat, delimiter='\t', quoting=csv.QUOTE_NONE)

	# master list of all dictionaries of positions of variants
	featcall_pos_list = []
	# unique set of encountered variants
	used_variants = set()

	### Loop through each chromosome / line of feat and creat a dictionary for each
	for line in feat_reader1:
		chr_dict = {}

		### in a Chr, loop through the variant calls for each gene/unit to gather their positions
		for index, call1 in enumerate(line[1].split(",")):
			variant = call1.split(":")[1]

			#skip variants if they are Unmapped or ^
			if variant == "Unmapped" or variant[-1] == "^":
				# print(f"{variant} not used", file = sys.stderr)
				continue

			#add passed variants to unique encountered set
			used_variants.add(variant)

			#if variant call is not in dictionary, create entry of a list of positions
			if variant not in chr_dict:
				chr_dict[variant] = [line[2].split(",")[index]]
			#else add a new position to the list for that variant
			else:
				chr_dict[variant].append(line[2].split(",")[index])

		# add the finished dictionary for this Chr to the master list
		featcall_pos_list.append(chr_dict)

# sort the used variants for looping
u_variants = list(used_variants)
u_variants.sort()

#####  Loop through used variants and measure how many Chr they are on and how localized they are within Chr
#####  outputting those that pass
### Loop through used variants
for u_variant in u_variants:
	# reset counters
	chrom_count = 0
	high_hits = 0
	mid_hits = 0
	low_hits = 0
	max_dist = 0
	
	### Loop through each chromosome
	for chrom in featcall_pos_list:
		chrom_count += 1

		if u_variant in chrom:
			### categorize occurences for each Chr based on how many 
			occurences = len(chrom[u_variant])

			print(f"{u_variant} is in chrom_{chrom_count} {occurences} times", file = sys.stderr)
			if occurences == 1:
				low_hits += 1
			elif occurences > 3:
				high_hits += 1
			else:
				mid_hits += 1

			### also calculate max distance between occurences on that Chr
			if occurences > 1:
				#intentionally leaves out first position as it is always 0 and not useful for distance
				#there is no good way to deal with the effectively unknown position of first entry (is limitation of current .feat)
				cur_dist = int(chrom[u_variant][-1]) - int(chrom[u_variant][1])
				if cur_dist > max_dist:
					max_dist = cur_dist

	### for each variant, output it if it's positioning passes filters
	if (low_hits == 0 and (mid_hits + high_hits) == 1) or (low_hits == 1 and mid_hits == 0 and high_hits == 1):
		if max_dist < distance_cutoff:
			print(u_variant)
		else:
			# print distance cutoffs info in case it is informative for choosing cutoff
			print(f"{u_variant} is {max_dist} bp, cutoff was {distance_cutoff}", file = sys.stderr)
	