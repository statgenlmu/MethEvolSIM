# Load required libraries
library(parallel)
library(optparse)
library(devtools); load_all()

# Define command-line options
option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Design file, full path", metavar = "character"),
  make_option("--output-path", type = "character", default = NULL,
              help = "Full path for output files", metavar = "character"),
  make_option("--n-cores", type = "integer", default = NULL,
              help = "Number of cores", metavar = "integer"),
  make_option("--n-tips", type = "integer", default = NULL,
              help = "Number of tree tips", metavar = "integer"),
  make_option("--CFTP", type = "logical", default = FALSE, 
              metavar = "TRUE/FALSE", help = "Use CFTP algorithm at tree root [default: %default]"),
  make_option(c("--server-name"), type = "character", default = NULL, 
              help = "Server name (fuego, chile, australia, mauritius, rio)")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("input", "output-path", "n-cores", "n-tips", "CFTP", "server-name")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

# Check if server name is valid
server_names <- c("fuego", "chile", "australia", "mauritius", "rio")
if (!(opt[["server-name"]] %in% server_names)) {
  stop("Invalid server name. Choose from: fuego, chile, australia, mauritius, rio")
}

# Load the simulation design file
load(opt[["input"]])

# Define index subsets
all_indices <- 1:nrow(sampled_params)
subsets <- split(all_indices, cut(all_indices, breaks = length(server_names), labels = server_names))

# Get the corresponding subset
indices <- subsets[[opt[["server-name"]]]]

# Generate pad_n based on the number of digits in n_sim
# To save the files with padded numbers so that are later listed in order with list.files()
pad_n <- nchar(as.character(all_indices)) + 1

# Print the indices for the selected server
print(paste("Simulating data in", opt[["server-name"]], "for the following indices:\n"))
print(indices)

# Save time of start for the subset
total_time_start <- Sys.time()

# Log file for failed cases
log_file <- paste0("failed_", opt[["server-name"]], ".log")
if (file.exists(log_file)) file.remove(log_file)  # Clear previous log
# Vector for saving failed indices
failed_indices <- c()

# Define the function to run each simulation
run_simulation <- function(i) {
  
  # Set the padded index for the output file name
  padded_index <- formatC(i, width = pad_n, format = "d", flag = "0")
  
  # Set a seed
  set.seed(i)
  
  # Simulate data
  output <- simulate_evolData(infoStr = spatial_str, tree = scaled_trees[[i]], params = sampled_params[i,], CFTP = opt[["CFTP"]])
  
  # Extract methylation data at the tree tips
  data <- list()
  for (tip in 1:opt[["n-tips"]]) {
    data[[tip]] <- output$data[[1]][[tip]]$seq
  }
  
  # Save the simulated data 
  save(data, file = file.path(opt[["output-path"]], paste0("abc_dataSIM_", padded_index, ".RData")))
}

# Run the simulations
mclapply(indices, function(i) {
  tryCatch({
    print(paste("Running simulation number:", i))
    start_time <- Sys.time()
    run_simulation(i)  
    end_time <- Sys.time()
    print(paste("Simulation", i, "finished in", difftime(end_time, start_time, units = "mins"), "mins."))
  }, error = function(e) {
    print(paste("Simulation", i, "failed. Error message in:", log_file))
    log_entry <- sprintf("Simulation %d failed: %s\n", i, e$message)
    write(log_entry, file = log_file, append = TRUE)
    failed_indices <- c(failed_indices, i)
    return(NULL)  # Return NULL to fake non-error
    #return(list(sim_id = i, error = e$message))  # Return failure info
  })
}, mc.cores = opt[["n-cores"]])

print("Done. Failed indices:")
print(failed_indices)

# Save time of end for the subset
total_time_end <- Sys.time()

# Print total simulation time
print(paste("Total simulation time:", difftime(total_time_end, total_time_start, units = "days"), "days."))
