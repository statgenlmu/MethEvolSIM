

library(devtools)
load_all()
library(parallel)
library(optparse)

# Define the command-line options
option_list <- list(
  make_option(c("-o", "--output-dir"), type = "character", default = NULL,
              help = "Full path to the output directory", metavar = "character"),
  make_option(c("-f", "--design-file"), type = "character", default = NULL,
              help = "Full path to the simulation design file", metavar = "character"),
  make_option(c("-b", "--branch-length"), type = "integer", default = NULL,
              help = "length of branch for simulating each evol step", metavar = "integer"),
  make_option(c("-s", "--start"), type = "integer", default = NULL,
              help = "start step number for simulating evol", metavar = "integer"),
  make_option(c("-e", "--end"), type = "integer", default = NULL,
              help = "end step number for simulating evol", metavar = "integer"),
  make_option(c("-p", "--paramComb-n"), type = "integer", default = NULL,
              help = "Number of parameter combinations", metavar = "integer"),
  make_option(c("-t", "--test-n"), type = "integer", default = NULL,
              help = "test number ", metavar = "integer"),
  make_option(c("-n", "--name-pattern"), type = "character", default = NULL,
              help = "output .RData files start name pattern", metavar = "character"),
  make_option(c("-r", "--replicate-n"), type = "integer", default = NULL,
              help = "number of data replicates to simulate", metavar = "integer")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("output-dir", "design-file", "branch-length", "start", "end", "paramComb-n", "test-n", "name-pattern", "replicate-n")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}


simul_CFTP_branch <- function(custom_params, index_params, b_length, start, end, out_digit_n, spatial_str, test_n, replicate_n){
  # Set the name for the output file with the padded parameter index
  padded_index_params <- formatC(index_params, width = 2, format = "d", flag = "0")
  padded_replicate_n <- formatC(replicate_n, width = 2, format = "d", flag = "0")
  out_file <- paste0(opt[["output-dir"]], "/", opt[["name-pattern"]], test_n, "_paramsID_", padded_index_params, "_rep_", replicate_n, ".out")
  # Redirect both the stout and stderr to the same file
  sink(out_file, type = c("output", "message"), append = TRUE)
  print(paste("Running CFTP_testConvergence: ", test_n, ". paramsID:", padded_index_params, ". Replicate:", padded_replicate_n))
  print("Given customized parameter values:")
  print(custom_params)
  if(start == 1){
    print("Generating combiStructureGenerator instance")
    combi <- combiStructureGenerator$new(infoStr = spatial_str, params = custom_params)
    # Save initial instance state and methylation data
    data <- list()
    for (str in 1:combi$get_singleStr_number()){
      data[[str]]<- transform_methStateEncoding(combi$get_singleStr(str)$get_seq())
    }
    padded_sim_n <- formatC(0, width = out_digit_n, format = "d", flag = "0")
    padded_replicate_n <- formatC(replicate_n, width = 2, format = "d", flag = "0")
    RData_file <- paste0(opt[["output-dir"]], "/", opt[["name-pattern"]], test_n, "_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_sim_n, ".RData")
    save(data, combi, file = RData_file)
    # Call cftp method from copy of initial instance, save instance state and methylation data
    print("Cloning and calling $cftp() method")
    cftp_combi <- combi$copy()
    cftp_combi$cftp()
    data <- list()
    for (str in 1:combi$get_singleStr_number()){
      data[[str]]<- transform_methStateEncoding(combi$get_singleStr(str)$get_seq())
    }
    RData_file <- paste0(opt[["output-dir"]], "/", opt[["name-pattern"]], test_n, "_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_cftp.RData")
    save(data, cftp_combi, file = RData_file)
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
    padded_replicate_n <- formatC(replicate_n, width = 2, format = "d", flag = "0")
    RData_file <- paste0(opt[["output-dir"]], "/", opt[["name-pattern"]], test_n, "_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_sim_n, ".RData")
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

# Set function to enable sampling initial eqFreqs according to given parameters
get_private <- function(x) {
  x[['.__enclos_env__']]$private
}

simul_CFTP_tests <- function(index_params){
  # Set seed
  set.seed(index_params)
  # Load the combinations of sampled parameters and spatial structure
  load(opt[["design-file"]])
  # Set parameter combination case
  custom_params <- sampled_params[index_params,]
  # Set IWE rate to 0
  custom_params$mu <- 0
  # Sample initial equilibrium frequencies for each structure
  samplereqFreqsI <- singleStructureGenerator$new("U", 1, params = custom_params)
  samplereqFreqsNI <- singleStructureGenerator$new("M", 1, params = custom_params)
  spatial_str$u_eqFreq <- rep(NA, dim(spatial_str)[1])
  spatial_str$p_eqFreq <- rep(NA, dim(spatial_str)[1])
  spatial_str$m_eqFreq <- rep(NA, dim(spatial_str)[1])
  for(str in 1:dim(spatial_str)[1]){
    if(spatial_str$globalState[str]=="U"){
      spatial_str[str,c("u_eqFreq", "p_eqFreq", "m_eqFreq")] <- get_private(samplereqFreqsI)$sample_eqFreqs()
    }
    
    if(spatial_str$globalState[str]=="M"){
      spatial_str[str,c("u_eqFreq", "p_eqFreq", "m_eqFreq")] <- get_private(samplereqFreqsNI)$sample_eqFreqs()
    }
  }
  # Call function to simulate data for each replicate
  for(r in 1:opt[["replicate-n"]]){
    simul_CFTP_branch(custom_params = custom_params,
                      index_params = index_params,
                      b_length = opt[["branch-length"]], # Set branch length
                      start = opt[["start"]], # Set start and end for the number of times to conduct simulations along the branch
                      end = opt[["end"]],
                      out_digit_n = 4, # Set number for maximum width of padded simulation number 
                      spatial_str = spatial_str,
                      test_n = opt[["test-n"]],
                      replicate_n = r)
  }
  
}

# Run in parallel using mclapply
mclapply(1:opt[["paramComb-n"]], function(index_params) simul_CFTP_tests(index_params), mc.cores = opt[["paramComb-n"]])



