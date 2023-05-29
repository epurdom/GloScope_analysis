This folder contains the code to run simulations with different parameters.

Details of simulation parameters please see the manuscript's Method section.

pre_process.R contains the code to create the reference data for simulation.

muscat_sim_serials.R contains the function we adapted from MUSCAT package to create simulations.

True_KL_serials_update.R contains the function to calculate the true distance for certain simulation settings.

Run ddir.R to create the folders for simulation results. 

To run certain simulation settings, please change the R filename in the example.sh to the corresponding R file name (e.g sd013_lfc005.R) in this folder, and submit the SLURM jobs.

After obtaining the simulation results, summary.R contains the code to create the csv files for power, ANOSIM test statistics.

boxplot.R contains the code to summarize the within or between group distance for boxplot visualization.


