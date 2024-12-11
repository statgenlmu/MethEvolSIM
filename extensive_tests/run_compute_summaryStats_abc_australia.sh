#!/bin/bash
#SBATCH --job-name=newSumStatsFreqPM_compute_summaryStats

# Compute summary statistics using the R script

Rscript /home/saracv/MasterThesis/MethEvolSIM/extensive_tests/compute_summaryStats.R \
  -d /scratch/saracv/abc/old_design \
  -p abc_dataSIM \
  -f /scratch/saracv/abc/old_design/abc_designSIM.RData \
  -n 100

