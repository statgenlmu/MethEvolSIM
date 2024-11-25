library(optparse)
# TODO: ADD when source is eliminated
# library(devtools); load_all()

## TODO: Once functions are commented and tested and belong to the package, this should be library(MethEvolSIM) / while developing library(devtools); load_all()
source("../../functions_summaryStats.R") # Load functions to compute summary statistics

if (FALSE){
  # Define the command-line options
  option_list <- list(
    make_option(c("-d", "--data-dir"), type = "character", default = NULL,
                help = "Full path to the directory containing the .RData files", metavar = "character"),
    make_option(c("-f", "--design-file"), type = "character", default = NULL,
                help = "Full path to the simulation design file", metavar = "character"),
    make_option(c("-s", "--stats"), type = "character", default = "all",
                help = "Comma-separated list of summary statistics to compute (default: all). Options: meanFreqP_i, meanFreqP_ni, sdFreqP_i, sdFreqP_ni, meanFreqM_i, meanFreqM_ni, sdFreqM_i, sdFreqM_ni, meanFracMoverMU_i, meanFracMoverMU_ni, sdFracMoverMU_i, sdFracMoverMU_ni, FChangeCherry_i, FChangeCherry_ni, Fitch, Transitions, meanCor", metavar = "character"),
    make_option(c("-n", "--sample-n"), type = "integer", default = NULL,
                help = "Number of samples per file", metavar = "integer"),
    make_option(c("-p", "--pattern"), type = "character", default = NULL,
                help = ".RData files start name pattern", metavar = "character"),
    make_option(c("-u", "--update-file"), type = "character", default = NULL,
                help = "Full path to existing summary statistics file", metavar = "character")
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

pad_n <- nchar(as.character(opt[["n-sim"]])) + 1

RData_files <- list.files(pattern = paste0("^abc_dataSIM_\\d{", pad_n, "}\\.RData"))

# Print the input files
print(RData_files)