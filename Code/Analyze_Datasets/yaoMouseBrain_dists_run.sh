#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

#SBATCH --mail-user=wtorous@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
R CMD BATCH --no-save Code/Analyze_Datasets/yaoMouseBrain_dists_run.R Code/Analyze_Datasets/yaoMouseBrain_dists_run.out
