library(optparse)
# TODO: ADD when source is eliminated
# library(devtools); load_all()

## TODO: Once functions are commented and tested and belong to the package, this should be library(MethEvolSIM) / while developing library(devtools); load_all()
source("../../functions_summaryStats.R") # Load functions to compute summary statistics

if (FALSE){
  # Define the command-line options
  option_list <- list(
    make_option("--data-dir", type = "character", default = NULL,
                help = "Full path to the directory containing the .RData files", metavar = "character"),
    make_option("--input-file", type = "character", default = NULL,
                help = "Name of the input .RData file", metavar = "character"),
    make_option("--design-file", type = "character", default = NULL,
                help = "Name of the simulation design file", metavar = "character"),
    make_option("--spatialStr-file", type = "character", default = NULL,
                help = "Name of the file with the spatial structure of the genomic region", metavar = "character"),
    make_option("--sample-n", type = "integer", default = NULL,
                help = "Number of samples per file", metavar = "integer"),
    make_option("--n-sim", type = "character", default = NULL,
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
                help = "Number of CpGs at each side of an island to exclude when computing mean correlations", metavar = "integer")
  )
}



# Define command-line options
option_list <- list(
  make_option("--n-sim", type = "integer", default = NULL,
              help = "Number of simulations", metavar = "integer")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
##TODO: Finish this

# Arguments as variables
dir <- opt[["data-dir"]]
input_file <- opt[["input-file"]]
sample_n <- opt[["sample-n"]]
u_threshold <- opt[["u-threshold"]]
m_threshold <- opt[["m-threshold"]]
cherry <- opt[["cherry-index"]]
minN_CpG <- opt[["minN-CpG"]]
shore_length <- opt[["shore-length"]]
pad_n <- nchar(as.character(length(opt[["n-sim"]]))) + 1
load(file.path(dir, opt[["design-file"]]))
load(file.path(dir, opt[["spatialStr-file"]]))

# Structural indeces of simulated CpG islands and non-islands
index_islands <- which(spatial_str$globalState=="U")
index_nonislands <- which(spatial_str$globalState=="M")


######################## Per sim function ######################################


#Function to compute summary statistics
compute_sumStats <- function(input_file, dir) {
  load(file.path(dir, input_file))  # Load the RData file (assumes `data` is stored in it)
  index <- as.integer(sub("abc_dataSIM_(\\d+)\\.RData", "\\1", input_file))
  
  # Compute example summary statistics
  MeanSiteFChange <- MeanSiteFChange_cherry(data, tree = scaled_trees[[index]], index_islands, index_nonislands)
  # Initialize the data frame to store the summary statistics
  sumStats <- data.frame(islandMeanFreqP = get_islandMeanFreqP(index_islands, data, sample_n),
                         islandSDFreqP = get_islandSDFreqP(index_islands, data, sample_n),
                         nonislandMeanFreqP = get_nonislandMeanFreqP(index_nonislands, data, sample_n),
                         nonislandSDFreqP = get_nonislandSDFreqP(index_nonislands, data, sample_n),
                         islandMeanFreqM = get_islandMeanFreqM(index_islands, data, sample_n),
                         islandSDFreqM = get_islandSDFreqM(index_islands, data, sample_n),
                         nonislandMeanFreqM = get_nonislandMeanFreqM(index_nonislands, data, sample_n),
                         nonislandSDFreqM = get_nonislandSDFreqM(index_nonislands, data, sample_n),
                         meanCor_i = compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 5, data, sample_n),
                         meanCor_ni = compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 5, data, sample_n),
                         MeanSiteFChange_i = MeanSiteFChange$island_meanFChange[cherry],
                         MeanSiteFChange_ni = MeanSiteFChange$nonisland_meanFChange[cherry],
                         cherryDist = MeanSiteFChange$dist[cherry],
                         Fitch_islandGlbSt = sum(computeFitch_islandGlbSt(index_islands, data, tree = scaled_trees[[index]], u_threshold, m_threshold)))
  
  # Define output filename
  
  # Set the padded index for the output file name
  padded_index <- formatC(index, width = pad_n, format = "d", flag = "0")
  output_file <- file.path(dir, paste0("sumStats_", padded_index, ".RData"))
  
  # Save summary statistics
  save(sumStats, file = output_file)
  
  # Print progress
  cat(sprintf("Processed: %s -> %s\n", input_file, output_file))
}

