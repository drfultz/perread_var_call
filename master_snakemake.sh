#/bin/bash
#PBS -N sm_master
#PBS -m abe
#PBS -M drfultz@iu.edu
#PBS -l nodes=1:ppn=1,vmem=2gb,walltime=48:00:00
#PBS -d /N/slate/drfultz/snakemake_test

TARGET=polished/Sha_teloB_edits_polished1X.fa

snakemake --cluster "qsub -m abe -M drfultz@iu.edu -e logs/{wildcards.sample}_{rule}.err -o logs/{wildcards.sample}_{rule}.out -l nodes=1:ppn={params.ppn},vmem={params.vmem},walltime={params.walltime}" --jobs 100 --use-conda $TARGET
