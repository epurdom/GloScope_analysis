#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1

#SBATCH --mail-user=wtorous@berkeley.edu
#SBATCH --mail-type=ALL

R CMD BATCH --no-save Code/Analyze_Datasets/pelkaCRC_dists_run.R Code/Analyze_Datasets/pelkaCRC_dists_run.out
