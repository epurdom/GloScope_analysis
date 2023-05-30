#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

#SBATCH --mail-user=wtorous@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
R CMD BATCH --no-save Code/Analyze_Datasets/melmsCOVIDLung_tester_run.R Code/Analyze_Datasets/melmsCOVIDLung_tester_run.out
