#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

#SBATCH --mail-user=wtorous@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
R CMD BATCH --no-save Code/Process_Datasets/process_ledergorMyeloma_default_scvi.R Code/Process_Datasets/process_ledergorMyeloma_default_scvi.out 
