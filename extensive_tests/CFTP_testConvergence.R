library(devtools)
load_all()
library(parallel)



simul_CFTP_branch <- function(custom_params, index_params, b_length, start, end, out_digit_n, spatial_str){
  # Set the name for the output file with the padded parameter index
  padded_index_params <- formatC(index_params, width = 2, format = "d", flag = "0")
  out_file <- paste0("/scratch/saracv/CFTP_test/CFTP_testConvergence_paramsID_", padded_index_params, ".out")
  # Redirect both the stout and stderr to the same file
  sink(out_file, type = c("output", "message"), append = TRUE)
  print(paste("Running CFTP_testConvergence_paramsID", padded_index_params))
  print("Given customized parameter values:")
  print(custom_params)
  if(start == 1){
    print("Generating combiStructureGenerator instance")
    combi <- combiStructureGenerator$new(infoStr = spatial_str, params = custom_params)
  }
  # Simulate evolution along branch n times
  print(paste("Simulating evolution along branch of length", b_length, end - start + 1, "times."))
  for (i in start:end){
    print(paste("Running simulation number:", i))
    # Simulate data
    combi$branch_evol(branch_length = b_length, dt = 0.01)
    # Save methylation data
    data <- list()
    for (str in 1:combi$get_singleStr_number()){
      data[[str]]<- transform_methStateEncoding(combi$get_singleStr(str)$get_seq())
    }
    # Save simulated data
    padded_sim_n <- formatC(i, width = out_digit_n, format = "d", flag = "0")
    if (i == end){
      # The last time, save also the combiStructureGenerator instance, to be able to start new simulations from last state
      save(data, combi, file = paste0("/scratch/saracv/CFTP_test/CFTP_testConvergence_paramsID_", padded_index_params, "_n_", padded_sim_n, ".RData" ))
    } else {
      save(data, file = paste0("/scratch/saracv/CFTP_test/CFTP_testConvergence_paramsID_", padded_index_params, "_", padded_sim_n, ".RData" ))
    }
  }
  # Stop redirecting output and messages
  sink()
}

simul_CFTP_tests <- function(index_params){
  # Set branch length
  b_length <- 10
  # Set start and end for the number of times to conduct simulations along the branch
  start <- 1
  end <- 50
  out_digit_n <- 4
  # Load the combinations of sampled parameters and spatial structure
  load("/scratch/saracv/abc_designSIM.RData")
  # Set parameter combination case
  custom_params <- sampled_params[index_params,]
  # Set IWE rate to 0
  custom_params$mu <- 0
  # Call function to simulate data
  simul_CFTP_branch(custom_params = custom_params,
                    index_params = index_params,
                    b_length = b_length,
                    start = start,
                    end = end,
                    out_digit_n = out_digit_n,
                    spatial_str = spatial_str)
  
}

# Run in parallel using mclapply
paramComb_n <- 10
mclapply(1:paramComb_n, function(index_params) simul_CFTP_tests(index_params), mc.cores = 10)




