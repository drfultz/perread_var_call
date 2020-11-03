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


paf_file = snakemake.input[0]
input_fq = snakemake.input[1]
min_align_length = snakemake.params.min_len
output_fq = snakemake.output[0]

#print(paf_file + readsbam + output_fq)

with open(paf_file,'r') as paf:
	paf_reader = csv.reader(paf, delimiter='\t', quoting=csv.QUOTE_NONE)
	aligned_names = set()
	for align in paf_reader:
		aligned_names.add(align[0])

with open(output_fq, "w") as outf:
	for aligned_name in aligned_names:
		print(aligned_name)
		pre_positions = [0]
		num_F_mappings = 0
		num_R_mappings = 0

		with open(paf_file,'r') as loop_paf:
			loop_paf_reader = csv.reader(loop_paf, delimiter='\t', quoting=csv.QUOTE_NONE)
			for aln in loop_paf_reader:
				if aligned_name == aln[0] and (int(aln[8]) - int(aln[7])) > min_align_length:
					if aln[4] == "-":
						num_R_mappings += 1
						pre_positions.append(int(aln[3]))
					elif aln[4] == "+":
						num_F_mappings += 1
						pre_positions.append(int(aln[2]))
					else:
						print("ERROR in paf")
						#this should probably exit at this point or have better handling

		positions = list(set(pre_positions))
		#print(positions)

		for read in pysam.FastxFile(input_fq):
			#print("\t" + read.name)
			count = 0
			#read_name = read.name
			if aligned_name == read.name.rsplit()[0]:
				read_length = len(read.sequence)
				positions.append(read_length)
				positions.sort()
				print(positions)
				#print(read.quality)

				orientation = ""
				oriented_sequence = ""
				oriented_qualities = ""
				oriented_positions = []

				if num_R_mappings > num_F_mappings:
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
				print(oriented_positions)

				for begin,end in zip(oriented_positions[:-1],oriented_positions[1:]):
					count += 1
					print(begin, "\t: ", end, "\tlength", len(oriented_sequence[begin:end]))

					a = pysam.FastxRecord()
					a.name = "%s:%s" % (aligned_name, str(count).zfill(2))
					a.sequence = oriented_sequence[begin:end]
					a.quality = oriented_qualities[begin:end]
					#print(a)
					outf.write(str(a) + "\n")
				print()
				break #assumes no duplicate readnames. Readname occurences after first will be skipped.
		else: # no break
			print("!!!!!!!!!!!!!!!!!!!!\nno match for {} found\n\n".format(aligned_name))
