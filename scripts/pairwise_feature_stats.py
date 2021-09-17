#!/usr/bin/env python3
import sys
import re
import csv

#tabulate pair counts for all provided features

features_file = sys.argv[1]
#feature_list_file = sys.argv[2]
table_file = "%s.feat" % features_file.rsplit('.', 1)[0]

# features_file = snakemake.input[0]
# feature_list_file = snakemake.input[1]
# table_file = snakemake.output[0]

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

table = Vividict() #assumes Python 3.7+ for ordering of dictionary's output
keys_set = set()

# with open(feature_list_file,'r') as feature_list:
# 	for entry in feature_list:
# 		keys_set.add(entry[:-1])
# 	#print(keys_set)
# 	keys_list = list(keys_set)
# 	keys_list.sort()
# 	for keyA in keys_list:
# 		table[keyA] = Vividict()
# 		for keyB in keys_list:
# 			table[keyA][keyB] = 0
# 		#print(list(table[keyA].values()))

with open(features_file,'r') as feature_calls:
	for line in feature_calls:
		features_field = line[:-1].split("\t")[1]
		features_field = features_field.replace("^","")
		units = features_field.split(",")
		#print(units)

		unit_lists = []
		for unit in units:
			unit_lists.append(unit.split("&"))

		for unitA,unitB in zip(unit_lists[:-1],unit_lists[1:]):
			#print(unitA,unitB)
			for featA in unitA:
				for featB in unitB:
					#print(featA," to ",featB)
					try:
						table[featA][featB] += 1
					except:
						table[featA][featB] = 1


#print(table)
print("feat/feat\t" + "\t".join([str(x) for x in table]))
for keyA in table:
	line = [keyA]
	for keyB in table:
		if table[keyA][keyB] == {}:
			line.append(0)
		else:
			line.append(table[keyA][keyB])
	print("\t".join([str(x) for x in line]))

# with open(table_file,"w") as out_table:
