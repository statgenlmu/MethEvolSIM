#!/bin/bash
#SBATCH --job-name=CFTP_testConvergence_local

# First test run: b_length set to 1, start to 1, end to 250.
# Second test run: b_length set to 0.1, start to 1, end to 500. 

# Define the test number variable
test_n=2

# Define the base directory path
base_dir="/scratch/saracv/CFTP_test/test${test_n}"

  Rscript CFTP_testConvergence.R -o $base_dir\
                                 -f /scratch/saracv/abc_designSIM.RData \
                                 -b 1 \
                                 -s 1 \
                                 -e 250 \
                                 -p 10 \
                                 -t $test_n \
                                 -n CFTP_testConvergence \
                                 -r 10
                                 
                                 


