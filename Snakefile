configfile: "config.yaml"
localrules: fa_to_fq, to_split_fq, call_splitfq, merge_calls, format_calls #,map_for_splitting

###minimap2 and paf-mappings_to-split-fq_sm.py can take .fa or .fq inputs and produce outputs which function downstream and maintain q scores, if present, in split fastqs.
###Using a duplicate for the mapping_for_split and splitter rules allows either .fa or .fq as inputs while still requiring a proper file extesion of either .fa or .fq

#temp rule
rule fa_to_fq:
    input:
        "input_sequence/{sample}.fa"
    output:
        "input_sequence/{sample}.fq"
    params:
        ppn="1",
        vmem="8gb",
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
        ppn="6",
        vmem="24gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c -t {params.ppn} -x map-ont {input[0]} {input[1]} > {output[0]}"

rule to_split_fq:
    input:
        "map_to_split/{sample}_to_split-ref.paf",
        "input_sequence/{sample}.fq"
    output:
        "oriented_sequence/{sample}_oriented.fq",
        "split_sequence/{sample}_split.fq"
    params:
        min_len=config["MIN_ALIGN_LENGTH"],
        ppn="1",
        vmem="12gb",
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
        ppn="2",
        vmem="16gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c --paf-no-hit -t {params.ppn} -x map-ont {input[1]} {input[0]} > {output[0]} \n"

rule call_splitfq:
    input:
        "split_mapping/{sample}_split_to_{feature}.cig_paf"
    output:
        "split_calls/{sample}_callsfor_{feature}.tsv"
    params:
        ppn="1",
        vmem="8gb",
        walltime="8:00:00"
    script:
        "scripts/multi_alignment_paf_caller_sm.py"

rule merge_calls:
    input:
        ["split_calls/{{sample}}_callsfor_{feat}.tsv".format(feat=f) for f in config["FEATURES"]]
    output:
        "merged_calls/{sample}_mergedcalls.tsv"
    params:
        ppn="1",
        vmem="8gb",
        walltime="1:00:00"
    shell:
        "paste {input} > {output}"

rule format_calls:
    input:
        "merged_calls/{sample}_mergedcalls.tsv"
    output:
        "feature_calls/{sample}.feat"
    params:
        ppn="1",
        vmem="8gb",
        walltime="1:00:00"
    script:
        "scripts/format_feature_calls_sm.py"

### The following rules are for quantifying transcripts or other whole reads
### (quantitative sequencing - qs) along a target .feat file
rule map_qs_to_feature:
    input:
        "features/{feature}.fa",
        "quant_seq/{qs_run}.fq"
    output:
        "qs_mapped/{qs_run}_to_{feature}.cig_paf"
    params:
        ppn="2",
        vmem="10gb",
        walltime="8:00:00"
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -c --secondary=no -t {params.ppn} -x map-ont {input[0]} {input[1]} > {output[0]}"

rule make_qs_track:
    input:
        "feature_calls/{sample}.feat",
        "qs_mapped/{qs_run}_to_{feature}.cig_paf",
        "qs_select_variants/{sel_var}_variants.txt"
    output:
        "qs_tracks/{sample}_with_{qs_run}_to_{feature}_at_{sel_var}.tsv"
    params:
        ppn="1",
        vmem="10gb",
        walltime="8:00:00"
    script:
        "scripts/qs_to_track_sm.py"

