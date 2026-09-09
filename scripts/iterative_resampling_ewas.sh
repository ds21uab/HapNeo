#!/bin/bash
## ============================================================================
## 05_iterative_resampling_gwas.sh
## SLURM array job: runs 05_iterative_resampling_gwas.R for iterations 1-100
## Submit with: sbatch 05_iterative_resampling_gwas.sh
## ============================================================================

#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH --partition=medium
#SBATCH --array=1-100
#SBATCH -e slurm-%A_%a.err
#SBATCH -o slurm-%A_%a.out
#SBATCH --job-name=hapneo_iter

date
echo "Running iteration ${SLURM_ARRAY_TASK_ID}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Create output directories
mkdir -p results/resampling


module load R/4.4.1-gfbf-2023b

Rscript iterative_resampling_ewas.R

date
