library(optparse)
library(ape)

# Define command-line options with unique letters
option_list <- list(
  make_option("--dir", type = "character", default = NULL,
              help = "Full path to save the output file", metavar = "character"),
  make_option("--file-name", type = "character", default = NULL,
              help = "Name of the output file", metavar = "character"),
  make_option("--CpG-n", type = "integer", default = NULL,
              help = "Number of CpG sites per genomic structure", metavar = "integer"),
  make_option("--str-n", type = "integer", default = NULL,
              help = "Number of structures of each type", metavar = "integer")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set the spatial structure
spatial_str <- data.frame(n = rep(opt[["CpG-n"]], 2*opt[["str-n"]]),
                          globalState = rep(c("U", "M"), opt[["str-n"]]))

# Save the spatial structure
file_path <- file.path(opt[["dir"]], opt[["file-name"]])
save(spatial_str, file = file_path)