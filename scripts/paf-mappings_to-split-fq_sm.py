#!/usr/bin/env python3
import pysam
import sys
import subprocess
import csv

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
min_align_length = snakemake.params.min_len
output_oriented_file = snakemake.output[0]
output_tagged_file = snakemake.output[1]

def generate_outputs():
	###Loop through the input fastq, looking for reads that match the alignment name
	for read in input_fq:
		if current_aln_name == read.name.split()[0]:
			print(current_aln_name)

			##Calculate and add the final position to positions
			read_length = len(read.sequence)
			pre_positions.append(read_length)

			##ensure that positions is only unique elements and sort
			positions = list(set(pre_positions)) 
			positions.sort()
			#print(positions)

			##Initialize
			orientation = ""
			oriented_sequence = ""
			oriented_qualities = ""
			oriented_positions = []

			##Reverse complement sequence and qualitites if needed
			if num_R_mappings > num_F_mappings: #if mostly R direction, reverse complement
				orientation = "RC"
				oriented_sequence = read.sequence.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
				if read.quality is not None:
					oriented_qualities = read.quality[::-1]
				oriented_positions = [read_length - x for x in positions][::-1]
			else:
				orientation = "unaltered"
				oriented_sequence = read.sequence
				if read.quality is not None:
					oriented_qualities = read.quality #"".join([chr(x + 33) for x in read.quality])
				oriented_positions = positions
			#print(orientation)
			#print(oriented_positions)

			##Write the whole original read in it's oriented state
			x = pysam.FastxRecord()
			x.name = current_aln_name
			x.sequence = oriented_sequence
			x.quality = oriented_qualities
			# print(x)
			output_oriented.write(str(x) + "\n")

			##Write the split reads (in the oriented direction) by looping through each pair of positions
			count = 0
			for begin,end in zip(oriented_positions[:-1],oriented_positions[1:]):
				count += 1
				#print(begin, "\t: ", end, "\tlength", len(oriented_sequence[begin:end]))
				a = pysam.FastxRecord()
				a.name = "%s:%s:%s" % (current_aln_name,str(count).zfill(2),begin)
				a.sequence = oriented_sequence[begin:end]
				a.quality = oriented_qualities[begin:end]
				#print(a)
				output_tagged.write(str(a) + "\n")

			##Readname occurences after first will be skipped.
			break

	###If no matching readname found while searching fastx
	else: # no break was found
		print(f"\n!!!!!!!!!!!!!!!!!!!!\nno match for {current_aln_name} found\n")


with open(paf_file,"r") as paf, open(output_oriented_file, "w") as output_oriented, open(output_tagged_file, "w") as output_tagged:
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	
	input_fq = list(pysam.FastxFile(input_fq_file)) # load fastq into RAM since it will be looped many times

	#####Initialize the record values
	current_aln_name = "ignore_this_message"
	pre_positions = [0]
	num_F_mappings = 0
	num_R_mappings = 0

	#####Loop through the paf once.
	for aln in paf_reader:

		#####If the alignment name has changed from the previous alignment record...
		if aln[0] != current_aln_name:

			###...generate the outputs for the previous group of alignments and write them...
			generate_outputs()

			###...then reset for the next group of alignments.
			current_aln_name = aln[0]
			pre_positions = [0]
			num_F_mappings = 0
			num_R_mappings = 0

		#####For the current alignment, if it is above the minimum alignment length...
		if (int(aln[8]) - int(aln[7])) > min_align_length:
			###...add to either reverse mappings...
			if aln[4] == "-":
				num_R_mappings += 1
				pre_positions.append(int(aln[3])) #for R mappings, use 'end' position
			###...or forward mappings.
			elif aln[4] == "+":
				num_F_mappings += 1
				pre_positions.append(int(aln[2])) #for F mappings use 'start' position
			else:
				raise Exception(f"\n\tImproperly formatted paf.\n\tField 5 should only contain '+' or '-' , but '{aln[4]}' was found")

	###Generate and write the output of the last group of alignments
	generate_outputs()
	print()