configfile: "config.yaml"
localrules: fa_to_fq, to_split_fq, call_splitfq, format_calls
wildcard_constraints: feature="|".join(config["FEATURES"])

###fasta input files are converted to fastq so that either format can go into the pipeline
rule fa_to_fq:
    input:
        "input_sequence/{sample}.fa"
    output:
        "input_sequence/{sample}.fq"
    params:
        ppn="1",
        memperCPU="8gb",
        walltime="1:00:00"
    conda:
        "envs/bbmap2.yaml"
    shell:
        "reformat.sh in={input[0]} out={output[0]}"

rule map_for_splitting:
    input:
        config["Splitting_Seq"],
        "input_sequence/{sample}.fq"
    output:
        "map_to_split/{sample}_to_split-ref.paf"
    params:
        ppn="3",
        memperCPU="4gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c --secondary=no -t {params.ppn} -x map-ont {input[0]} {input[1]} > {output[0]}"

rule to_split_fq:
    input:
        "map_to_split/{sample}_to_split-ref.paf",
        "input_sequence/{sample}.fq"
    output:
        "oriented_sequence/{sample}_oriented.fq",
        "split_sequence/{sample}_split.fq"
    params:
        min_len=config["MIN_ALIGN_LENGTH"],
        reorient=config["REORIENT"],
        ppn="1",
        memperCPU="8gb",
        walltime="8:00:00"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/paf-mappings_to-split-fq_sm.py"

rule map_splitfq:
    input:
        "split_sequence/{sample}_split.fq",
        "features/{feature}.fa"
    output:
        "split_mapping/{sample}_split_to_{feature}.cig_paf"
    params:
        ppn="8",
        memperCPU="3gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c --paf-no-hit -t {params.ppn} -x map-ont {input[1]} {input[0]} > {output[0]}"

rule call_splitfq:
    input:
        "split_mapping/{sample}_split_to_{feature}.cig_paf"
    output:
        "split_calls/{sample}_callsfor_{feature}.tsv"
    params:
        VARIANT_LEN_COV_MIN=config["VARIANT_LEN_COV_MIN"],
        VARIANT_LEN_COV_BEST=config["VARIANT_LEN_COV_BEST"],
        ppn="1",
        memperCPU="8gb",
        walltime="8:00:00"
    script:
        "scripts/multi_alignment_paf_caller_sm.py"

rule format_calls:
    input:
        expand("split_calls/{{sample}}_callsfor_{feat}.tsv", feat=config["FEATURES"])
    output:
        "feature_calls/{sample}.feat"
    params:
        ppn="1",
        memperCPU="8gb",
        walltime="1:00:00"
    script:
        "scripts/format_feature_calls_sm.py"

### The following rules are for quantifying transcripts or other whole reads
### (qs stands for quantitative sequencing) along a target .feat file
rule map_qs_to_feature:
    input:
        "features/{feature}.fa",
        "quant_seq/{qs_run}.fq"
    output:
        "qs_mapped/{qs_run}_to_{feature}.cig_paf"
    params:
        ppn="4",
        memperCPU="3gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c --secondary=no -t {params.ppn} -x sr {input[0]} {input[1]} > {output[0]}"

rule make_qs_track:
    input:
        "feature_calls/{sample}.feat",
        "qs_mapped/{qs_run}_to_{feature}.cig_paf",
        "qs_select_variants/{sel_var}_variants.txt"
    output:
        "qs_tracks/{sample}_with_{qs_run}_to_{feature}_at_{sel_var}.tsv"
    params:
        ppn="1",
        memperCPU="10gb",
        walltime="8:00:00"
    script:
        "scripts/qs_to_track_sm.py"

