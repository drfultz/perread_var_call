## perread_var_call
Long arrays of tandem repeats such as 45S/48S, 5S, and pericentromeric sequences are common in genomes and are impossible to assemble using traditional methods - short reads cannot be used, and most long read assembly software collapses these regions.  This pipeline uses existing tools and custom scripts to examine long DNA reads from repetitive regions for variant or 'morph' content which can be used downstream to assemble.

The pipeline takes a fasta or fastq (.fa or .fq) of long reads as input.  For each read in the input, it outputs an entry in a feature calls file (.feat file).  Each feature calls entry contains information on the structural variants (indels, VNTR, etc. or presence/absence) present within each repeat unit, in order, within a read.  Thus, each DNA sequence read is translated into its higher order sequence of variants.  These sequences of variants can then be used to assemble repeat regions using the sister pipeline `varcall_longread_assembler` (https://github.com/drfultz/varcall_longread_assembler).  This pipeline was built and used to resolve the 45S rDNA arrays in _Arabidopsis thaliana_ from ultra-long ONT reads, but was built to take custom inputs for general use in tandem repeats.  Alternatively, this pipeline can be used for simple analysis of the ordering of repeat variants in contigs or reference sequence from repetitive regions.  A few scripts are included in the `tools` folder for simple analysis of sequences of repeat variants from .feat files.

An earlier version of this pipeline was used in the first near complete assembly of the _A. thaliana_ NOR/rDNA regions https://www.science.org/doi/full/10.1126/sciadv.adj4509 and this sequence is now included in the reference genome (Col-CC/TAIR12) https://phoenixbioinfo.org/announcing-the-publication-of-tair12-consult-the-fully-reannotated-genome-on-the-european-nucleotide-archive/

### Getting started
The pipeline is designed to be run by `snakemake` using its built-in interface with the package manager `conda`.  If needed, install conda, then use conda to install `snakemake`.  `snakemake` will install all other dependencies when the pipeline is first run with the `--use-conda` option (as below).  Basic knowledge of `snakemake` usage will be helpful in working with this pipeline.  This pipeline was tested with `snakemake` v7.32.4
```
conda create -n snakemake -c conda-forge -c bioconda snakemake=7.32.4
```

Once you have `snakemake`, install a copy of the project:
```
git clone https://github.com/drfultz/perread_var_call.git
cd perread_var_call
snakemake -j 1 --use-conda feature_calls/reads-test.feat
```
Included in the repository are test inputs.  It will likely take several minutes for conda to install all of the packages.  The pipeline itself should only take a minute or two to run.  Upon completion, `feature_calls/reads-test.feat` should be written and contain 13 lines of feature calls (which correspond to a subset of the input .fq entries that contained repeat sequence).

## Running your own variant calling
### Basic setup
1. Copy your reads (`.fa` or `.fq` extension only) into the `input_sequence/` folder, for example `input_sequence/SAMPLE.fq`.
2. For each region of interest within a repeat, place one `.fa` file (containing one entry for each known variant) in the `features/` folder.  For ease of use later, keep the names simple, for example `A.fa` instead of  `3ETS_regionA.fa`.
3. Place a reference sequence into the `refs/` folder.  This will be used for dividing repeats units.
4. Update the `config.yaml`.
    - A) Under the `FEATURES:` entry, enter one line for each file in the `features/` folder that you want to call.  Use the filename without the extension.  For example, enter the line `    - A` to tell the pipeline to use `features/A.fa` to call variants.
    - B) After `Splitting_Seq:` enter the path to your reference from step 3. For example `refs/MYREF.fa`.
    - C) Adjust the integer after `MIN_ALIGN_LENGTH:`. This value should be smaller than the length of the Splitting_Seq.
5. From the base directory, run the pipeline:
	```
	snakemake -j 4 --use-conda feature_calls/SAMPLE.feat
	```

### Pipeline Outputs
The key outputs of the pipeline are the `.feat` file present in the `feature_calls/` folder and the `_oriented.fq` file present in the `oriented_sequence/` folder.  For both outputs, only reads with matches to the `Splitting_Seq` will be kept.
#### Feature calls (`.feat` file) format
Each line of a .feat file represents the feature calls for a retained read from the input sequence file.  Feature calls are either variant calls or presence vs absence calls, depending on the contents of the .fa file for that 'feature' in the `features/` folder.  Each line is divided into three fields separated by tabs: 1. readname, 2. feature calls, 3. unit start sites.
1) The readname field is the associated readname from the input .fa or .fq.
2) The feature calls field is divided into individual repeat units by a comma, with units in the order of the oriented read (see below).  Within each unit, there are calls for each user-defined feature which are joined by '&'.  A call is made for each feature, as ordered in the `config.yaml`, regardless of the order within the DNA sequence of the unit.  Calls for each feature are formatted with the feature name, a colon, then the call, for example `A:VAR1` is a VAR1 variant at the region defined as A.  For example, here are the calls for a three unit long read for features 'A' and 'B':
	```
	A:VAR1&B:458,A:VAR2&B:294,A:VAR1%&B:Unmapped
	```
	Read segments that only partially cover a variant or have no mapping are given the call 'Unmapped'. The cutoff mapping length for this can be adjusted in the `config.yaml` `VARIANT_LEN_COV_MIN`. Mappings that cover the mapped variant at a higher threshold have the call followed by a '%'. The cutoff mapping length for this can be adjusted in the `config.yaml` `VARIANT_LEN_COV_BEST`.
3) The unit start sites field represents the start positions of each repeat unit (as defined by `Splitting_Seq`) for the oriented read.  Commas separate each successive start site.

#### `_oriented.fq` file
The oriented fastq file contains reads that were variant called (any read that has the `Splitting_Seq`).  The reads (and feature calls in the .feat) are reoriented in the forward direction relative to the majority of the mappings to the `Splitting_Seq`, unless the `REORIENT` field in the `config.yaml` is changed.  (Reads from regions with repeat units facing both directions are unlikely to be oriented consistently, as this pipeline was designed for direct repeats.)

## Optimizing variant calling
#### Input
While this pipeline currently does not utilize fastq quality values, it will maintain them in the `_oriented.fq` output.  Otherwise, there is no benefit to either fasta or fastq input.  For efficient performance, it is strongly recommended to first filter the input to only the desired reads, for example by presence of the target repeats and by read length.

#### Suggestions for features sequence design
For each feature, each variant sequence should be flanked by sequence shared between all variants.  For example, fasta entries for a 100bp insertion would include 1) a 'reference'  sequence without the insertion, consisting of the 75bp of sequence on either side of the point of the insertion and 2) an insertion variant sequence with the same set of 75bp sequences on either side of the 100bp insertion.  This ensures that both variant sequences contain the same shared sequence, which leaves the presence or absence of the insert as the determining factor for which is chosen as the variant call.  This also applies to tandem repeats or more complex repeat structures, where the same invariable sequence should flank each variant.

Alternatively, calls can be made for presence/absence of a sequence by using only one sequence in a `.fa` in the `features/` folder but this may be slightly less accurate.

For both applications keep in mind that the default settings of the mapping program `minimap2` and its alignment score will be used to map to the variants and thus determine the variant calls.

#### Choice of Splitting_Seq and MIN_ALIGN_LENGTH
The `Splitting_Seq` .fa file is used as the point where reads are divided into repeat units (the first nucleotide of each alignment to the `Splitting_Seq` in a read is the start point for each segment).  Thus, it is important for this sequence to be well conserved in all repeats and free of variants, so that it can be recognized consistently in low accuracy reads.  The `MIN_ALIGN_LENGTH` parameter should probably be in the hundreds of bp or higher, depending on repeat size, repeat conservation, and the users desire for sensitivity vs. accuracy when splitting.  The sequence in the `Splitting_Seq` .fa file can be much longer than the `MIN_ALIGN_LENGTH` if desired.  This would cause any match within this region to be used to split even if the rest of the region is absent.  Multiple entries in the `Splitting_Seq` .fa can also be used, unit splits will be made at a mapping to any entry.

### Experimental
Code for quantifying variants from quantitative sequencing (RNA-seq, DiMeLo, etc.) is included but is not well tested. (The last two rules in the Snakefile.)

### Limitations
Currently, this pipeline does not classify the direction of individual units in reads.  For repeat regions with alternating directions of repeats, this pipeline will therefore have limited value.  Similarly, the structure within individual units (for example the ordering of each of the examined features/SVs) is not considered as it is assumed that different units are all structured in the same fashion.  The ability to recognize the direction and structure of repeats will be added in a v2 of this pipeline, which is being actively developed.

## Pipeline Diagram
![dag](./dag.svg)