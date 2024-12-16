library(optparse)
library(parallel)
## TODO: Modify how the package is loaded before submitting
library(devtools); load_all()

# Define command-line options
option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Design file, full path", metavar = "character"),
  make_option("--n-cores", type = "integer", default = NULL,
              help = "Number of cores", metavar = "integer"),
  make_option("--n-tips", type = "integer", default = NULL,
              help = "Number of tree tips", metavar = "integer"),
  make_option("--CFTP", type = "logical", default = FALSE, 
    metavar = "TRUE/FALSE", help = "Use CFTP algorithm at tree root [default: %default]")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load the simulation design file
load(opt[["input"]])
print("here")
print(head(sampled_params))
print(scaled_trees)
print(spatial_str)

# Define the number of runs and cores
n_sim <- nrow(sampled_params)
n_cores <- opt[["n-cores"]]
print("here 2")
print(opt[["n-cores"]])
print(opt[["CFTP"]])

# Generate pad_n based on the number of digits in n_sim
# To save the files with padded numbers so that are later listed in order with list.files()
pad_n <- nchar(as.character(n_sim)) + 1

run_simulation_subset <- function(start_idx, end_idx) {
  
  # Generate the padded start and end indexes for writing the output
  padded_start <- formatC(start_idx, width = pad_n, format = "d", flag = "0")
  padded_end <- formatC(end_idx, width = pad_n, format = "d", flag = "0")
  
  # Set the name for the .out file for the entire subset
  out_file <- paste0("abc_dataSIM_", padded_start, "_to_", padded_end, ".out")
  
  # Open the log file once for the entire subset
  log_connection <- file(out_file, open = "w")
  
  # Redirect both regular output and message output to the same file
  sink(log_connection, type = "output")
  sink(log_connection, type = "message")
  
  # Run the simulations for the current subset
  for (sim in start_idx:end_idx) {
    print(paste("Running simulation number:", sim))
    
    padded_index <- formatC(sim, width = pad_n, format = "d", flag = "0")
    
    # Set a seed
    set.seed(sim)
    
    # Simulate data
    output <- simulate_evolData(infoStr = spatial_str, tree = scaled_trees[[sim]], params = sampled_params[sim,], CFTP = opt[["CFTP"]])
    
    # Extract methylation data at the tree tips
    data <- list()
    for (tip in 1:opt[["n-tips"]]) {
      data[[tip]] <- output$data[[1]][[tip]]$seq
    }
    # Save the simulated data
    save(data, file = paste0("abc_dataSIM_", padded_index, ".RData"))
    
  }
  
  # Stop redirecting output and messages
  sink(type = "output")
  sink(type = "message")
  
  # Close the connection to the log file
  close(log_connection)
}


print("here3")

# Calculate the number of runs per chore
runs_per_chore <- ceiling(n_sim / n_cores)

# Divide the runs among the cores
run_indices <- split(seq_len(n_sim), ceiling(seq_along(seq_len(n_sim)) / runs_per_chore))

# Create a list of start and end indices for each chunk
chunks <- lapply(run_indices, function(idx) c(min(idx), max(idx)))

# Run the simulation subsets in parallel
mclapply(chunks, function(chunk) run_simulation_subset(chunk[1], chunk[2]), mc.cores = n_cores)
