#!/bin/bash
#SBATCH --job-name=newSumStatsFreqPM_compute_summaryStats

# Compute summary statistics. Script requires R library optparse
# Currently needs to be run from git repo MethEvolSIM branch develop, under directory extensive tests or similar. 

Rscript /home/saracv/MasterThesis/MethEvolSIM/extensive_tests/compute_summaryStats.R -d /scratch/saracv \ # Full path to the directory containing the .RData files
										     -p abc_dataSIM_0 \ # Full path to the simulation design file
                                                                                     -f /scratch/saracv/abc_designSIM.RData \  # .RData files start name pattern
                                                                                     -n 100 # Number of samples per file (samples = tips = replicates)
                                                                                     #-s # Comma-separated list of summary statistics to compute (default: all). Options: meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanFracMoverMU_i,meanFracMoverMU_ni,sdFracMoverMU_i,sdFracMoverMU_ni,FChangeCherry_i,FChangeCherry_ni,Fitch,Transitions,meanCor
                                                                                     #-u Full path to existing summary statistics file (option to update selection of summary statistics)

