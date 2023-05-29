#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mail-user=hao_wang@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-100

R CMD BATCH --no-save sd013_lfc005.R simulation_$SLURM_ARRAY_TASK_ID


