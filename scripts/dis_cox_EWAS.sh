#!/bin/bash
#SBATCH -N 1                      # Number of nodes
#SBATCH --cpus-per-task=10       # Number of CPUs
#SBATCH --mem=100G               # Memory
#SBATCH -t 12:00:00              # Walltime
#SBATCH --partition=short       # Partition
#SBATCH -e slurm-%A_%a.err       # Error file per array task
#SBATCH -o slurm-%A_%a.out       # Output file per array task

date

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load R/4.4.1-gfbf-2023b

# Call your R script with array index as argument
Rscript dis_cox_EWAS.R $SLURM_ARRAY_TASK_ID

date
