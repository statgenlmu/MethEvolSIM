#!/bin/bash
#SBATCH --job-name=CFTP_compute_summaryStats

for i in {1..10}
do
  if [ $i -lt 10 ]
  then
    param_id="nestedReps_CFTP_testConvergence2_paramsID_0$i"
  else
    param_id="nestedReps_CFTP_testConvergence2_paramsID_$i"
  fi
    
    # Call the R script with the new rep_id
    Rscript compute_summaryStats.R -d /scratch/saracv/CFTP_test/test2 \
                                   -p $param_id \
                                   -f /scratch/saracv/abc_designSIM.RData \
                                   -n 10 \
                                   -s meanFreqP_i,meanFreqP_ni,sdFreqP_i,sdFreqP_ni,meanFreqM_i,meanFreqM_ni,sdFreqM_i,sdFreqM_ni,meanCor
  done
done



