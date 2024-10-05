#!/bin/bash
#SBATCH --job-name=CFTP_compute_summaryStats

# Define the test number variable
test_n=2

# Define the base directory path
base_dir="/scratch/saracv/CFTP_test/test${test_n}"

#for i in $(seq 1 10)
for i in 2 3 4 5 7 9 10
do
  if [ $i -lt 10 ]
  then
    param_id="nestedReps_CFTP_testConvergence${test_n}_paramsID_0$i"
  else
    param_id="nestedReps_CFTP_testConvergence${test_n}_paramsID_$i"
  fi
    
    # Print param_id to check its value
    echo "param_id is set to: $param_id"
  
    # Call the R script for $branch_evol output data
    Rscript compute_summaryStats.R -d $base_dir \
                                   -p $param_id \
                                   -f /scratch/saracv/abc_designSIM.RData \
                                   -n 10 \
                                   -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Steepness
                                   
    # Call the R script for $cftp output data
    Rscript compute_summaryStats_cftpOutput.R -d $base_dir \
    					      -p $param_id \
    					      -f /scratch/saracv/abc_designSIM.RData \
    					      -n 10 \
    					      -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Steepness	
done



