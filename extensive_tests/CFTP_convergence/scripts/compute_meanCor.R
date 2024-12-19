## TODO: modify how the package is loaded before submitting
library(devtools)
load_all()
library(parallel)
library(optparse)

# Define the command-line options
option_list <- list(
  make_option("--dir", type = "character", default = NULL,
              help = "Full path to dir for input and output files", metavar = "character"),
  make_option("--design-file", type = "character", default = NULL,
              help = "Full path to the simulation design file", metavar = "character"),
  make_option("--n-sim", type = "integer", default = NULL,
              help = "Number of simulations", metavar = "integer"),
  make_option("--replicate-n", type = "integer", default = NULL,
              help = "number of data replicates to simulate", metavar = "integer"),
  make_option(c("-s", "--start"), type = "integer", default = NULL,
              help = "start step number for simulating evol", metavar = "integer"),
  make_option(c("-e", "--end"), type = "integer", default = NULL,
              help = "end step number for simulating evol", metavar = "integer")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("dir", "n-sim", "replicate-n", "design-file", "start", "end")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

print(opt[["dir"]])
print(opt[["n-sim"]])
print(opt[["replicate-n"]])
print(opt[["design-file"]])
print(opt[["start"]])
print(opt[["end"]])


# Define function to calculate mean correlation for data simulated with branch steps
computeCor_branchSteps <- function(index_params, params_pad_n, start, end, replicate_n, replicate_pad_n, step_pad_n, dir){
    print(paste("Computing mean correlation at branch steps parameter ID:", index_params))
    padded_index_params <- formatC(index_params, width = params_pad_n, format = "d", flag = "0")
    for(rep in 1:replicate_n){
        print(paste("Replicate number:", rep))
        padded_replicate_n <- formatC(rep, width = replicate_pad_n, format = "d", flag = "0")
    
        # Initialize the data frame to store the summary statistics
        summaryStats <- data.frame(meanCor_i = rep(NA, length((start-1):end)), meanCor_ni = rep(NA, length((start-1):end)))
    
        for (step in (start-1):end){
            padded_step_n <- formatC(step, width = step_pad_n, format = "d", flag = "0")
            load(file.path(dir, paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_step_n, ".RData")))
            meanCor_i[step + 1] <- compute_meanCor_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data, sample_n = 1)
            meanCor_ni[step + 1] <- compute_meanCor_i(index_nonislands = index_nonislands, minN_CpG = 80, shore_length = 10, data, sample_n = 1)
        }

        # Save the data frame for each replicate
        save(summaryStats, file.path(dir, paste0("summaryStats_CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_step_n, ".RData")))
    }
}

# Define function to compute mean correlation for data simulated with the CFTP method
computeCor_CFTP <- function(index_params, params_pad_n, replicate_n, replicate_pad_n, dir){
    print(paste("Computing mean correlation after CFTP parameter ID:", index_params))
    padded_index_params <- formatC(index_params, width = params_pad_n, format = "d", flag = "0")

    # Initialize the data frame to store the summary statistics
    summaryStats <- data.frame(meanCor_i = rep(NA, replicate_n), meanCor_ni = rep(NA, replicate_n))
    for (rep in 1:replicate_n){
        print(paste("Replicate number:", rep))
        padded_replicate_n <- formatC(rep, width = replicate_pad_n, format = "d", flag = "0")
        load(file.path(dir, paste0("CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_cftp.RData")))
        meanCor_i[rep] <- compute_meanCor_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data, sample_n = 1)
        meanCor_ni[rep] <- compute_meanCor_i(index_nonislands = index_nonislands, minN_CpG = 80, shore_length = 10, data, sample_n = 1)
        print("here")
        
    }
    # Save the data frame for each parameter index
    save(summaryStats, file.path(dir, paste0("summaryStats_CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, "_", padded_step_n, "_cftp.RData")))
}

# Define function to compute meanCor for both branch steps and CFTP for a given parameter combination (index)
compute_meanCor <- function(index_params){
    padded_index_params <- formatC(index_params, width = params_pad_n, format = "d", flag = "0")
    out_file <- file.path(dir, paste0("compute_meanCor_CFTP_testConvergence_paramsID_", padded_index_params, ".out"))

    # Open the log file once for the given parameters and replicate number
    #log_connection <- file(out_file, open = "w")
  
    # Redirect both regular output and message output to the same file
    #sink(log_connection, type = "output")
    #sink(log_connection, type = "message")
    
    # Compute the mean correlations
    computeCor_branchSteps(index_params, params_pad_n, start, end, replicate_n, replicate_pad_n, step_pad_n, dir)
    computeCor_CFTP(index_params, params_pad_n, replicate_n, replicate_pad_n, dir)

    # Stop redirecting output and messages
    #sink()
}

# Load info simulation design: sampled_params, spatial_str, tree
load(opt[["design-file"]])

# Structural indeces of simulated CpG islands and non-islands
index_islands <- which(spatial_str$globalState=="U")
index_nonislands <- which(spatial_str$globalState=="M")

# Set the parameters for the functions
dir = opt[["dir"]]
params_pad_n <- nchar(as.character(opt[["n-sim"]]))
start = opt[["start"]]
end = opt[["end"]]
replicate_n = opt[["replicate-n"]]
replicate_pad_n = nchar(as.character(opt[["replicate-n"]])))
step_pad_n = nchar(as.character(opt[["end"]])) + 1

# Try out
load("/scratch/saracv/CFTP_testConvergence/SSE_imp/design.RData")
index_islands <- which(spatial_str$globalState=="U")
index_nonislands <- which(spatial_str$globalState=="M")
dir = "/scratch/saracv/CFTP_testConvergence/SSE_imp"
params_pad_n <- nchar(as.character(10))
start = 1
end = 251
replicate_n = 10
replicate_pad_n = nchar(as.character(5))
step_pad_n = nchar(as.character(2)) + 1
n_sim = 10

# Run in parallel using mclapply

compute_meanCor(index_params = 1)

mclapply(1:n_sim, function(index_params) compute_meanCor(index_params), mc.cores = n_sim)
