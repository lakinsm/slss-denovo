#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=cpu_compute
#SBATCH --cpus-per-task=48
#SBATCH --job-name=asffast_nextflow_sbatch
#SBATCH --mem=160gb
#SBATCH --time=10:00:00
#SBATCH --output=asffast_nextflow_sbatch_%j.log

module load python37
export NXF_CLUSTER_SEED=$( shuf -i 0-16777216 -n 1 )

echo -e "\n$(date): Enqueuing job"
srun
echo -e "\n$(date): Job finished for data directory ${2}\n"
