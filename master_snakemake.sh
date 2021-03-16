#/bin/bash
#PBS -N sm_master
#PBS -m abe
#PBS -M NAME@iu.edu
#PBS -l nodes=1:ppn=1,vmem=8gb,walltime=30:00:00
#PBS -d /N/slate/drfultz/working_perread_var_call

TARGET=formatted_calls/test.feat

source ~/.bashrc
conda activate snakemake

snakemake --cluster "qsub -m abe -M NAME@iu.edu -e logs/{wildcards.sample}_{rule}.err -o logs/{wildcards.sample}_{rule}.out -l nodes=1:ppn={params.ppn},vmem={params.vmem},walltime={params.walltime}" --jobs 100 --use-conda --local-cores 1 $TARGET
