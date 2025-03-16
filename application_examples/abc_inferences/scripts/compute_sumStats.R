library(optparse)
library(parallel)
# TODO: ADD when source is eliminated
library(devtools); load_all()

# Define the command-line options
option_list <- list(
  make_option("--data-dir", type = "character", default = NULL,
              help = "Full path to the directory containing the .RData files", metavar = "character"),
  make_option("--design-file", type = "character", default = NULL,
              help = "Name of the simulation design file", metavar = "character"),
  make_option("--spatialStr-file", type = "character", default = NULL,
              help = "Name of the file with the spatial structure of the genomic region", metavar = "character"),
  make_option("--sample-n", type = "integer", default = NULL,
              help = "Number of samples per file", metavar = "integer"),
  make_option("--n-sim", type = "integer", default = NULL,
              help = "Total number of simulations", metavar = "integer"),
  make_option("--u-threshold", type = "numeric", default = NULL,
              help = "Threshold for category unmethylated", metavar = "numeric"),
  make_option("--m-threshold", type = "numeric", default = NULL,
              help = "Threshold for category methylated", metavar = "numeric"),
  make_option("--cherry-index", type = "integer", default = NULL,
              help = "Index for cherry to use", metavar = "integer"),
  make_option("--minN-CpG", type = "integer", default = NULL,
              help = "Minimum number of CpGs to compute mean correlations", metavar = "integer"),
  make_option("--shore-length", type = "integer", default = NULL,
              help = "Number of CpGs at each side of an island to exclude when computing mean correlations", metavar = "integer"),
  make_option("--n-cores", type = "integer", default = NULL,
              help = "Number of cores", metavar = "integer"),
  make_option("--categorized-data", type = "logical", default = FALSE,
              help = "Methylation states are categorized", metavar = "logical")
  )



# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("data-dir", "design-file", "spatialStr-file", "sample-n", "n-sim", "u-threshold", "m-threshold", "cherry-index", "minN-CpG", "shore-length", "categorized-data")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

# Arguments as variables
dir <- opt[["data-dir"]]
sample_n <- opt[["sample-n"]]
u_threshold <- opt[["u-threshold"]]
m_threshold <- opt[["m-threshold"]]
cherry <- opt[["cherry-index"]]
minN_CpG <- opt[["minN-CpG"]]
shore_length <- opt[["shore-length"]]
categorized_data <- opt[["categorized-data"]]
pad_n <- nchar(as.character(opt[["n-sim"]])) + 1

load(file.path(dir, opt[["design-file"]]))
load(file.path(dir, opt[["spatialStr-file"]]))

# Structural indeces of simulated CpG islands and non-islands
index_islands <- which(spatial_str$globalState=="U")
index_nonislands <- which(spatial_str$globalState=="M")

files <- list.files(dir, pattern = "^abc_dataSIM_.*\\.RData$")

######################## Per sim function ######################################

#Function to compute summary statistics
compute_sumStats <- function(input_file, dir) {
  load(file.path(dir, input_file))  # Load the RData file 
  index <- as.integer(sub("abc_dataSIM_(\\d+)\\.RData", "\\1", input_file))
  
  # Compute example summary statistics
  MeanSiteFChange <- MeanSiteFChange_cherry(data = data, categorized_data = categorized_data, tree = scaled_trees[[index]], index_islands = index_islands, index_nonislands = index_nonislands)
  # Initialize the data frame to store the summary statistics
  sumStats <- data.frame(islandMeanFreqP = get_islandMeanFreqP(index_islands = index_islands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         islandSDFreqP = get_islandSDFreqP(index_islands = index_islands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         nonislandMeanFreqP = get_nonislandMeanFreqP(index_nonislands = index_nonislands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         nonislandSDFreqP = get_nonislandSDFreqP(index_nonislands = index_nonislands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         islandMeanFreqM = get_islandMeanFreqM(index_islands = index_islands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         islandSDFreqM = get_islandSDFreqM(index_islands = index_islands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         nonislandMeanFreqM = get_nonislandMeanFreqM(index_nonislands = index_nonislands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         nonislandSDFreqM = get_nonislandSDFreqM(index_nonislands = index_nonislands, data = data, categorized_data = categorized_data, sample_n = sample_n),
                         meanCor_i = compute_meanCor_i(index_islands = index_islands, minN_CpG = minN_CpG, shore_length = shore_length, data = data, sample_n = sample_n, categorized_data = categorized_data),
                         meanCor_ni = compute_meanCor_ni(index_nonislands = index_nonislands, minN_CpG = minN_CpG, shore_length = shore_length, data = data, sample_n = sample_n, categorized_data = categorized_data),
                         MeanSiteFChange_i = MeanSiteFChange$island_meanFChange[cherry],
                         MeanSiteFChange_ni = MeanSiteFChange$nonisland_meanFChange[cherry],
                         cherryFreqsChange_i = mean_CherryFreqsChange_i(data = data, categorized_data = categorized_data, index_islands = index_islands, tree = scaled_trees[[index]], pValue_threshold = 0.05)[1,"FreqsChange"],
                         treeFreqsChange_i = mean_TreeFreqsChange_i(tree = scaled_trees[[index]], data = data, categorized_data = categorized_data, index_islands = index_islands, pValue_threshold = 0.05)
  )
  
  # Define output filename
  
  # Set the padded index for the output file name
  padded_index <- formatC(index, width = pad_n, format = "d", flag = "0")
  output_file <- file.path(dir, paste0("sumStats_", padded_index, ".RData"))
  
  # Save summary statistics
  save(sumStats, file = output_file)
  
  # Print progress
  cat(sprintf("Processed: %s -> %s\n", input_file, output_file))
}


# Run the simulations
mclapply(files, function(i) {
  compute_sumStats(i, dir)
}, mc.cores = opt[["n-cores"]])



