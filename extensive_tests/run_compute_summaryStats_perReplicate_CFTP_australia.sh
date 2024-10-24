#!/bin/bash
#SBATCH --job-name=CFTP_compute_summaryStats

# Define the test number variable
test_n=1

# Define the base directory path
base_dir="/scratch/saracv/CFTP_test/test${test_n}"

# Outer loop for paramsID
for i in 1 2 3 4 5 6 7 8 9 10
do
  if [ $i -lt 10 ]
  then
    param_id_base="CFTP_testConvergence${test_n}_paramsID_0$i"
  else
    param_id_base="CFTP_testConvergence${test_n}_paramsID_$i"
  fi

  # Inner loop for replicates
  for r in 1 2 3 4 5 6 7 8 9 10
  do
    if [ $r -lt 10 ]
    then
      param_id="${param_id_base}_rep_0$r"
    else
      param_id="${param_id_base}_rep_$r"
    fi
    
    # Print param_id to check its value
    echo "param_id is set to: $param_id"
    
    # Call the R script for $branch_evol output data
    Rscript compute_summaryStats.R -d $base_dir \
                                   -p $param_id \
                                   -f /scratch/saracv/abc_designSIM.RData \
                                   -n 1 \
                                   -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Transitions
                                   
    # Call the R script for $cftp output data
    Rscript compute_summaryStats_cftpOutput.R -d $base_dir \
                                   -p $param_id \
                                   -f /scratch/saracv/abc_designSIM.RData \
                                   -n 1 \
                                   -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Transitions

  done  # Closing inner loop
done  # Closing outer loop




