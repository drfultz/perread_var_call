#!/usr/bin/env python3
import pysam
import sys
import subprocess
import csv


### Version for var_call snakemake.

### Takes a paf alignment file and uses the beginning position of target mappings to split
### fastq reads of matching names into subreads.

### Inputs are expected to be 1) a .paf file from a minimap2 alignment of the reads (as query)
### to the 'splitting sequence' (as reference) and 2) an unmapped bam of the reads to split.

### The output keeps only reads present in the alignment file.  The 'splitting sequence' is included
### in the outputs such that the output split fastq entries for a given input read
### add up to the full length of the input read.
### Reads 

### Intended to split long multi-repeat ONT reads into single rDNA units by first aligning them to a 25S reference


# paf_file = sys.argv[1]
# input_fq = sys.argv[2]
# min_align_length = 1000 #alignment length may need to be altered for each reference
# output_oriented_file = input_fq.split('.')[0] + "_orient.fq"
# output_tagged_file = input_fq.split('.')[0] + "_split.fq"

paf_file = snakemake.input[0]
input_fq = snakemake.input[1]
min_align_length = snakemake.params.min_len
output_oriented_file = snakemake.output[0]
output_tagged_file = snakemake.output[1]

#print(paf_file + input_fq + output_tagged_file)

with open(paf_file,"r") as paf, open(output_oriented_file, "w") as output_oriented, open(output_tagged_file, "w") as output_tagged:
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	#####Initialize the record values
	current_readname = "ignore_this_message"
	pre_positions = [0]
	num_F_mappings = 0
	num_R_mappings = 0

	for aln in paf_reader:
		if aln[0] != current_readname:
			#####Look for a match in the fastq to write as outputs of the old record
			for read in pysam.FastxFile(input_fq):
				#print("\t" + read.name)
				count = 0
				if current_readname == read.name.split()[0]:
					###Add final position (read length) to positions
					print(current_readname)
					read_length = len(read.sequence)
					pre_positions.append(read_length)
					positions = list(set(pre_positions)) #reduce positions to unique elements
					positions.sort()
					print(positions)

					###Initial values
					orientation = ""
					oriented_sequence = ""
					oriented_qualities = ""
					oriented_positions = []

					###Reverse complement sequence and qualitites if needed
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
					print(orientation)
					#print(oriented_positions)

					###Write the whole original read in it's oriented state
					x = pysam.FastxRecord()
					x.name = current_readname
					x.sequence = oriented_sequence
					x.quality = oriented_qualities
					# print(x)
					output_oriented.write(str(x) + "\n")

					###Write the split reads in the oriented direction
					for begin,end in zip(oriented_positions[:-1],oriented_positions[1:]):
						count += 1
						print(begin, "\t: ", end, "\tlength", len(oriented_sequence[begin:end]))

						a = pysam.FastxRecord()
						a.name = "%s:%s:%s" % (current_readname,str(count).zfill(2),begin)
						a.sequence = oriented_sequence[begin:end]
						a.quality = oriented_qualities[begin:end]
						#print(a)
						output_tagged.write(str(a) + "\n")
					print()

					###Readname occurences after first will be skipped.
					break #assumes no duplicate readnames in fastx

			###If no matching readname found while searching fastx
			else: # no break was found
				print(f"!!!!!!!!!!!!!!!!!!!!\nno match for {current_readname} found\n\n")

			#####Reset the record
			current_readname = aln[0]
			pre_positions = [0]
			num_F_mappings = 0
			num_R_mappings = 0


		if (int(aln[8]) - int(aln[7])) > min_align_length:
			if aln[4] == "-":
				num_R_mappings += 1
				pre_positions.append(int(aln[3]))
			elif aln[4] == "+":
				num_F_mappings += 1
				pre_positions.append(int(aln[2]))
			else:
				raise Exception(f"\n\tImproperly formatted paf.\n\tField 5 should only contain '+' or '-' , but '{aln[4]}' was found")
				#this should probably exit at this point or have better handling