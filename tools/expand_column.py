#!/usr/bin/env python3
import sys

### Takes a file with a single column of values (such as a single column
### of variant calls from a .feat file) and expands that into a grid
### with each column representing presence/absence of each variant found in the original.

if len(sys.argv)!=2:
	print("\n\tERROR: please input a single column of values.\n", file=sys.stderr)
	sys.exit(2)
	
column_file = sys.argv[1]

### read in file
with open(column_file, 'r') as column_obj:
	column = [line.rstrip() for line in column_obj]

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