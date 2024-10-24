#!/bin/bash
#SBATCH --job-name=newSumStatsFreqPM_compute_summaryStats

# Compute again the subset of summary statistics for which the functions have been updated

Rscript /home/sara/methylation/MethEvolSIM/extensive_tests/compute_summaryStats.R -d /home/sara/methylation/phd_project_saracv/simulation_studies/parallel_abc \
										   -p abc_dataSIM_0 \
                                                                                   -f /home/sara/methylation/phd_project_saracv/simulation_studies/abc_designSIM.RData \
                                                                                   -n 100 \
                                                                                   -s meanFreqM_i,meanFreqM_ni \
                                                                                   -u /home/sara/methylation/MethEvolSIM/extensive_tests/summaryStats_abc_dataSIM_0.RData


