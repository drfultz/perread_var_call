#!/usr/bin/env python3
import sys

### Takes a .feat file and a specified feature, and expands each found variant of the specified feature
### into a table where each row represents a unit and each column represents presence/absence of
### the specified variant at that unit.
### Multiple .feat file lines are merged into each other, and not differentiated

if len(sys.argv)!=3:
	print("\n\tERROR: please input a .feat file and a name of the feature.\n", file=sys.stderr)
	sys.exit(2)
	
feat_file = sys.argv[1]
target_feature = sys.argv[2]

### read in .feat and collect the calls for the chosen feature, one per unit
column = []
with open(feat_file, 'r') as feat_obj:
	for line in feat_obj:
		# skip empty lines	
		if not line.strip():
			continue

		### parse the units into a list
		units = line.rstrip().split("\t")[1].split(",")

		### loop the units and calls and extract the desired feature calls
		for unit in units:
			for call in unit.split("&"):
				name, variant = call.split(":", 1)
				if name == target_feature:
					column.append(variant)

### If no matches found, error
if not column:
	print(f"\n\tERROR: no calls found for feature '{target_feature}'.\n", file=sys.stderr)
	sys.exit(1)

### Reduce to set and sort
fields_set = list(set(column))
fields_set.sort()

###make header
print("\t".join(fields_set))

for entry in column:
	out_line = []
	for u in fields_set:
		if u == entry:
			out_line.append("1")
		else:
			out_line.append("0")

	print("\t".join(out_line))