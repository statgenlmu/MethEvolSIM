#!/bin/bash
#SBATCH --job-name=newSumStatsFreqPM_compute_summaryStats

# Compute again the subset of summary statistics for which the functions have been updated

Rscript /home/saracv/MasterThesis/MethEvolSIM/extensive_tests/compute_summaryStats.R -d /scratch/saracv \
										     -p abc_dataSIM_0 \
                                                                                     -f /scratch/saracv/abc_designSIM.RData \
                                                                                     -n 100 \
                                                                                     -s meanFreqM_i,meanFreqM_ni \
                                                                                     -u /home/saracv/MasterThesis/MethEvolSIM/extensive_tests/summaryStats_abc_dataSIM_0.RData

