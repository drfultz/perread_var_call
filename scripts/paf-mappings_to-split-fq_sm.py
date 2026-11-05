#!/usr/bin/env python3
import pysam
import csv
import sys
from collections import defaultdict

### Takes a paf alignment file and uses the beginning position of target mappings to split
### fastq reads of matching names into subreads.

### Inputs are expected to be 1) a .paf file from a minimap2 alignment of the reads (as query)
### to the 'splitting sequence' (as reference) and 2) a fastq of the reads to split.

### The output keeps only reads present in the alignment file.  The 'splitting sequence' is included
### in the outputs such that the output split fastq entries for a given input read
### add up to the full length of the input read.

### Intended to split long multi-repeat ONT reads into single rDNA units by first aligning them to a 25S reference

# paf_file = sys.argv[1]
# input_fq_file = sys.argv[2]
# min_align_length = 1000 #alignment length may need to be altered for each reference
# output_oriented_file = input_fq_file.split('.')[0] + "_orient.fq"
# output_tagged_file = input_fq_file.split('.')[0] + "_split.fq"

paf_file = snakemake.input[0]
input_fq_file = snakemake.input[1]
output_oriented_file = snakemake.output[0]
output_tagged_file = snakemake.output[1]
min_aln_len = snakemake.params.min_len
reorient = snakemake.params.reorient

def write_outputs(read, hit, reorient, output_oriented, output_tagged):
	##Calculate and add the final position to positions
	read_length = len(read.sequence)

	##ensure that positions is only unique elements and sort
	positions = sorted(set(hit["positions"] + [read_length])) 

	##Initialize
	oriented_sequence = ""
	oriented_qualities = ""
	oriented_positions = []

	##Reverse complement sequence and qualitites if set to do so and there are more R mappings
	if reorient == "Y" and  hit["rev"] > hit["fwd"]: #if mostly R direction, reverse complement
		comment = "RC"
		oriented_sequence = read.sequence.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
		if read.quality is not None:
			oriented_qualities = read.quality[::-1]
		oriented_positions = [read_length - x for x in positions][::-1]
	else:
		comment = "F"
		oriented_sequence = read.sequence
		if read.quality is not None:
			oriented_qualities = read.quality
		oriented_positions = positions

	##Write the whole original read in it's oriented state
	orient_read = pysam.FastxRecord()
	orient_read.name = read.name
	orient_read.comment = comment
	orient_read.sequence = oriented_sequence
	orient_read.quality = oriented_qualities
	# print(orient_read)
	output_oriented.write(str(orient_read) + "\n")

	##Write the split reads (in the oriented direction) by looping through each pair of positions
	count = 0
	for begin,end in zip(oriented_positions[:-1],oriented_positions[1:]):
		count += 1
		read_seg = pysam.FastxRecord()
		read_seg.name = f"{read.name}:{count:0{width}d}:{begin}"
		read_seg.sequence = oriented_sequence[begin:end]
		read_seg.quality = oriented_qualities[begin:end]
		#print(read_seg)
		output_tagged.write(str(read_seg) + "\n")


with open(paf_file,"r") as paf, open(output_oriented_file, "w") as output_oriented, open(output_tagged_file, "w") as output_tagged:
	##### Read in paf and process into dictionary of start positions
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	hits_dict = defaultdict(lambda: {"positions": [0], "fwd": 0, "rev": 0})
	for aln in paf_reader:
		if int(aln[8]) - int(aln[7]) < min_aln_len:
			continue
		h = hits_dict[aln[0]]
		if aln[4] == "+":
			h["fwd"] += 1
			h["positions"].append(int(aln[2]))
		else:
			h["rev"] += 1
			h["positions"].append(int(aln[3]))

	### Test for the longest number of units in a read to set zeroes used
	max_units = max((len(h["positions"]) for h in hits_dict.values()), default=1)
	width = len(str(max_units))

	##### Look for readname matches in hits from the fastq, and output when found
	for read in pysam.FastxFile(input_fq_file):
		hit = hits_dict.pop(read.name, None)
		if hit is not None:
			write_outputs(read, hit, reorient, output_oriented, output_tagged)

	### Output readnames that exist in hits (above min_aln_len) but not in fastq
	for name in hits_dict:
		print(f"no match for {name} found", file=sys.stderr)