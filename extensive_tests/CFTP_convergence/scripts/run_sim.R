## TODO: modify how the package is loaded before submitting
library(devtools)
load_all()
library(parallel)
library(optparse)


# Define the command-line options
option_list <- list(
  make_option("--dir", type = "character", default = NULL,
              help = "Full path to save the output file", metavar = "character"),
  make_option(c("-f", "--design-file"), type = "character", default = NULL,
              help = "Full path to the simulation design file", metavar = "character"),
  make_option(c("-b", "--branch-length"), type = "numeric", default = NULL,
              help = "length of branch for simulating each evol step", metavar = "numeric"),
  make_option(c("-s", "--start"), type = "integer", default = NULL,
              help = "start step number for simulating evol", metavar = "integer"),
  make_option(c("-e", "--end"), type = "integer", default = NULL,
              help = "end step number for simulating evol", metavar = "integer"),
  make_option(c("-r", "--replicate-n"), type = "integer", default = NULL,
              help = "number of data replicates to simulate", metavar = "integer")
)



# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("dir", "design-file", "branch-length", "start", "end", "replicate-n")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}



load(opt[["design-file"]])

n_sim <- nrow(sampled_params)

# Generate pad_n based on the number of digits in the number or simulations and the number of replicates
# To save the files with padded numbers so that are later listed in order with list.files()
params_pad_n <- nchar(as.character(n_sim))




simul_CFTP_branch <- function(custom_params, params_pad_n, index_params, b_length, start, end, step_pad_n, spatial_str, replicate_n, replicate_pad_n){
  
  # Set the name for the output file with the padded parameter index
  padded_index_params <- formatC(index_params, width = params_pad_n, format = "d", flag = "0")
  padded_replicate_n <- formatC(replicate_n, width = replicate_pad_n, format = "d", flag = "0")
  out_file <- file.path(opt[["dir"]], paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, ".out"))

  # Open the log file once for the given parameters and replicate number
  log_connection <- file(out_file, open = "w")
  
  # Redirect both regular output and message output to the same file
  sink(log_connection, type = "output")
  sink(log_connection, type = "message")

  print(paste("Running CFTP_testConvergence. paramsID:", padded_index_params, ". Replicate:", padded_replicate_n))
  print("Given customized parameter values:")
  print(custom_params)
  
  # Set seed
  set.seed(index_params)
  
  if(start == 1){
    
    print("Generating combiStructureGenerator instance")
    combi <- combiStructureGenerator$new(infoStr = spatial_str, params = custom_params)
    
    ## Save initial instance state and methylation data ##
    data <- list()
    for (str in 1:combi$get_singleStr_number()){
      data[[str]]<- transform_methStateEncoding(combi$get_singleStr(str)$get_seq())
    }
    
    padded_step_n <- formatC(0, width = step_pad_n, format = "d", flag = "0")
    RData_file <- file.path(opt[["dir"]], paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_step_n, ".RData"))
    save(data, combi, file = RData_file)
    
    ## Call cftp method from copy of initial instance, save instance state and methylation data ##
    print("Calling $cftp() method")
    cftp_combi <- combi$cftp()
    
    data <- list()
    for (str in 1:cftp_combi$get_singleStr_number()){
      data[[str]]<- transform_methStateEncoding(cftp_combi$get_singleStr(str)$get_seq())
    }
    
    RData_file <- file.path(opt[["dir"]], paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_cftp.RData"))
    save(data, cftp_combi, file = RData_file)
  }
  
  ## Simulate evolution along branch n times ##
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
    padded_step_n <- formatC(i, width = step_pad_n, format = "d", flag = "0")
    RData_file <- file.path(opt[["dir"]], paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_step_n, ".RData"))
    
    if (i == end){
      # The last time, save also the combiStructureGenerator instance, to be able to start new simulations from last state
      save(data, combi, file = RData_file)
    } else {
      save(data, file = RData_file)
    }
  }
  
  # Stop redirecting output and messages
  sink()
}



simul_CFTP_tests <- function(index_params){
  
  # Set parameter combination case
  custom_params <- sampled_params[index_params,]
  
  # Call function to simulate data for each replicate
  for(r in 1:opt[["replicate-n"]]){
    simul_CFTP_branch(custom_params = custom_params,
                      index_params = index_params,
                      params_pad_n = params_pad_n,
                      b_length = opt[["branch-length"]], # Set branch length
                      start = opt[["start"]], # Set start and end for the number of times to conduct simulations along the branch
                      end = opt[["end"]],
                      step_pad_n = nchar(as.character(opt[["end"]])) + 1, # Set number for maximum width of padded step number 
                      spatial_str = spatial_str,
                      replicate_n = r,
                      replicate_pad_n = nchar(as.character(opt[["replicate-n"]])))
  }
  
}

# Run in parallel using mclapply

mclapply(1:(opt[["end"]]- opt[["start"]]), function(index_params) simul_CFTP_tests(index_params), mc.cores = n_sim)






