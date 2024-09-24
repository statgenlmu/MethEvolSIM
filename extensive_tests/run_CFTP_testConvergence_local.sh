#!/bin/bash
#SBATCH --job-name=CFTP_testConvergence_local

# First test run: b_length set to 10, start to 1, end to 50
# Second test run: b_length set to 1, start to 1, end to 500

  Rscript CFTP_testConvergence.R -o /home/sara/methylation/MethEvolSIM/extensive_tests/test2 \
                                 -f /home/sara/methylation/phd_project_saracv/simulation_studies/abc_designSIM.RData \
                                 -b 1 \
                                 -s 1 \
                                 -e 3 \
                                 -p 10 \
                                 -t 2 \
                                 -n CFTP_testConvergence 
                                 
                                 


