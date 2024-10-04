#!/bin/bash
#SBATCH --job-name=CFTP_testConvergence_local

# First test run: b_length set to 10, start to 1, end to 50
# Second test run: b_length set to 0.1, start to 1, end to 500. 
# Debug run (test 3): commented out method $cftp. Leads to no changes
# Debut run (test 4): b_length set to 1, start to 1, end to 250
# Debug run (test 5): keeps previous, but uncomments $cftp
# Debug run (test 6): Modified -b in script from integer to numeric. b_length set to 0.1, start to 1, end to 500.

  Rscript CFTP_testConvergence.R -o /scratch/saracv/CFTP_test/test6 \
                                 -f /scratch/saracv/abc_designSIM.RData \
                                 -b 0.1 \
                                 -s 1 \
                                 -e 500 \
                                 -p 10 \
                                 -t 6 \
                                 -n CFTP_testConvergence \
                                 -r 10
                                 
                                 


