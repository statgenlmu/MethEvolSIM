#!/bin/bash
#SBATCH --job-name=newSumStatsFreqPM_compute_summaryStats

# Compute again the subset of summary statistics for which the functions have been updated

Rscript /home/saracv/MasterThesis/MethEvolSIM/extensive_tests/compute_summaryStats.R -d /scratch/saracv \
										     -p abc_dataSIM_0 \
                                                                                     -f /scratch/saracv/abc_designSIM.RData \
                                                                                     -n 100 \
                                                                                     -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor


