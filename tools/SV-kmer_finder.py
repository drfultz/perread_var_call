#!/usr/bin/env python3
import sys

#tabulate pair counts for all provided features

features_file = sys.argv[1]

selfpair_file = "%s_self-pair.tsv" % features_file.rsplit('.', 1)[0]
neighbor_file = "%s_neighbor.tsv" % features_file.rsplit('.', 1)[0]

### Define Vividict to allow new keys to be added on the fly
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

selfpair_vdict = Vividict() #assumes Python 3.7+ for ordering of dictionary's output
neighbor_vdict = Vividict()

with open(features_file,'r') as feature_calls:
	for line in feature_calls:
		features_field = line.strip().split("\t")[1]
		features_field = features_field.replace("^","")
		units = features_field.split(",")
		#print(units)

		unit_lists = []
		for unit in units:
			unit_lists.append(unit.split("&"))

		## Iterate within each unit for self-pair analysis
		for unitZ in unit_lists:
			for featX in unitZ:
				for featY in unitZ:
					try:
						selfpair_vdict[featX][featY] += 1
					except:
						selfpair_vdict[featX][featY] = 1			

		## Iterate each pair of units for neighbors analysis
		for unitA,unitB in zip(unit_lists[:-1],unit_lists[1:]):
			#print(unitA,unitB)
			for featA in unitA:
				for featB in unitB:
					#print(featA," to ",featB)
					try:
						neighbor_vdict[featA][featB] += 1
					except:
						neighbor_vdict[featA][featB] = 1


### Write the self-pairings file
with open(selfpair_file,"w") as selfpair_out:
	selfpair_out.write("feat/feat\t" + "\t".join([str(x) for x in selfpair_vdict]) + "\n")
	for keyA in selfpair_vdict:
		line = [keyA]
		for keyB in selfpair_vdict:
			if selfpair_vdict[keyA][keyB] == {}:
				line.append(0)
			else:
				line.append(selfpair_vdict[keyA][keyB])
		selfpair_out.write("\t".join([str(x) for x in line]) + "\n")

### Write the neighbors file
with open(neighbor_file,"w") as neighbor_out:
	neighbor_out.write("feat/feat\t" + "\t".join([str(x) for x in neighbor_vdict]) + "\n")
	for keyA in neighbor_vdict:
		line = [keyA]
		for keyB in neighbor_vdict:
			if neighbor_vdict[keyA][keyB] == {}:
				line.append(0)
			else:
				line.append(neighbor_vdict[keyA][keyB])
		neighbor_out.write("\t".join([str(x) for x in line]) + "\n")
