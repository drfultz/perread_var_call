#!/bin/bash
#SBATCH -J sm_master
#SBATCH -p general
#SBATCH -o sm_master_%j.out
#SBATCH -e sm_master_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drfultz@iu.edu
#SBATCH -D /N/slate/drfultz/working_perread_var_call
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=36:00:00

TARGET=formatted_calls/allColto2021-04-01_125kb+35-125kb-QC_map45S.feat

source ~/.bashrc
conda activate snakemake

snakemake --cluster 'sbatch \
-J sm_master \
-p general \
-o logs/{wildcards.sample}_{rule}_${{SLURM_JOBID}}.out \
-e logs/{wildcards.sample}_{rule}_${{SLURM_JOBID}}.err \
--mail-type=ALL \
--mail-user=drfultz@iu.edu \
--nodes=1 \
--ntasks=1 \
--cpus-per-task={params.ppn} \
--mem-per-cpu={params.vmem} \
--time={params.walltime} \
' --jobs 100 --use-conda --local-cores 1 $TARGET
