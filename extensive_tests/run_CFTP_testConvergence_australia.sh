#!/bin/bash
#SBATCH --job-name=CFTP_testConvergence_local

# First test run: b_length set to 0.1, start to 1, end to 500.  

# Define the test number variable
test_n=1

# Define the base directory path
base_dir="/scratch/saracv/CFTP_test/test${test_n}"

  Rscript CFTP_testConvergence.R -o $base_dir\
                                 -f /scratch/saracv/abc_designSIM.RData \
                                 -b 0.1 \
                                 -s 1 \
                                 -e 1000 \
                                 -p 10 \
                                 -t $test_n \
                                 -n CFTP_testConvergence \
                                 -r 10
                                 
                                 


