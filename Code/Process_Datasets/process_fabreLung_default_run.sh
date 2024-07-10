#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

#SBATCH --mail-user=hao_wang@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
R CMD BATCH --no-save Code/Process_Datasets/process_fabreLung_default_run.R Code/Process_Datasets/process_fabreLung_default.out
