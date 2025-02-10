library(optparse)

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input CpG_count.RData path", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output full path with file name", metavar = "character")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("input", "output")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

# Load CpG counts
load(opt[["input"]])

# Create the data frame spatial_str
spatial_str <- data.frame(n = CpG_count_region$CpG_count)

# Add the globalState column based on CpG_count_region$str values
spatial_str$globalState <- ifelse(CpG_count_region$str == "Island", "U", "M")

# Save spatial structure
save(spatial_str, file = opt[["output"]])