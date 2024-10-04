#!/bin/bash
#SBATCH --job-name=CFTP_compute_summaryStats

#for i in $(seq 1 10)
for i in 2 3 4 5 6 7 9 10
do
  if [ $i -lt 10 ]
  then
    param_id="nestedReps_CFTP_testConvergence2_paramsID_0$i"
  else
    param_id="nestedReps_CFTP_testConvergence2_paramsID_$i"
  fi
    
    # Call the R script for $branch_evol output data
    Rscript compute_summaryStats.R -d /scratch/saracv/CFTP_test/test5 \
                                   -p $param_id \
                                   -f /scratch/saracv/abc_designSIM.RData \
                                   -n 10 \
                                   -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Steepness
                                   
    # Call the R script for $cftp output data
    Rscript compute_summaryStats_cftpOutput.R -d /scratch/saracv/CFTP_test/test5 \
    					      -p $param_id \
    					      -f /scratch/saracv/abc_designSIM.RData \
    					      -n 10 \
    					      -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor,Steepness	
done



