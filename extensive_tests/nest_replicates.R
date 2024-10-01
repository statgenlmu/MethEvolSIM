

# Define the number of cases, reps, and steps
n_cases <- 10    # Number of cases (paramsID_01, paramsID_02, ...)
n_reps <- 10    # Number of replicates (e.g., rep1, rep2, ..., rep10)
n_steps <- 251  # Number of steps (from 0000 to 0250, inclusive)

# Loop over each case
for (case in 1:n_cases) {
  print(paste("Processing parameter combination case:", case))
  
  # Loop over each time step
  for (step in 0:(n_steps - 1)) {
    print(paste("Time step:", step))
    
    # Format the time step with leading zeros
    step_formatted <- sprintf("%04d", step)
    
    # Initialize an empty list to store the data for each case and time step
    rep_data <- list()
    
    # Loop over each replicate
    for (rep in 1:n_reps) {
      
      # Format the file name with the case, rep, and step
      file_name <- sprintf("/scratch/saracv/CFTP_test/test2/CFTP_testConvergence2_paramsID_%02d_rep_%02d_%s.RData", case, rep, step_formatted)
      
      # Load the .RData file
      load(file_name)
      
      # Append the data to rep_data list
      rep_data[[rep]] <- data  
    }
    out_file <- sprintf("/scratch/saracv/CFTP_test/test2/nestedReps_CFTP_testConvergence2_paramsID_%02d_%s.RData", case, step_formatted)
    data <- rep_data
    save(data, file = out_file)
  }
  
  # Initialize an empty list to store the data for each case 
  rep_data <- list()
  
  # Loop over $cftp() output files
  for(rep in 1:n_reps){
    # Format the file name with the case, rep, and step
    file_name <- sprintf("/scratch/saracv/CFTP_test/test2/CFTP_testConvergence2_paramsID_%02d_rep_%02d_cftp.RData", case, rep)
    
    # Load the .RData file
    load(file_name)
    
    # Append the data to rep_data list
    rep_data[[rep]] <- data 
  }
  out_file <- sprintf("/scratch/saracv/CFTP_test/test2/nestedReps_CFTP_testConvergence2_paramsID_%02d_cftp.RData", case)
  data <- rep_data
  save(data, file = out_file)
}

