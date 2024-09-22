#!/bin/bash
#SBATCH --job-name=CFTP_compute_summaryStats

Rscript compute_summaryStats.R -d /home/sara/methylation/MethEvolSIM/extensive_tests -f /home/sara/methylation/phd_project_saracv/simulation_studies/abc_designSIM.RData -n 1 -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor
