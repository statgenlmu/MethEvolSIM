library(optparse)

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input CpG_count.RData path", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output full path with file name", metavar = "character"),
  make_option(c("--name-selected-region"), type = "character", default = NULL,
              help = "Name of selected genomic region", metavar = "character"),
  make_option(c("--chr"), type = "character", default = NULL,
              help = "Chromosome ID of selected region", metavar = "character"),
  make_option(c("--start"), type = "integer", default = NULL,
              help = "Start index of selected genomic region", metavar = "integer"),
  make_option(c("--end"), type = "integer", default = NULL,
              help = "End index of selected genomic region", metavar = "integer"),
  make_option(c("--side-bp"), type = "integer", default = NULL,
              help = "Number of bp to select structures starting after start - side-bp and ending before end + side-bp", metavar = "integer")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("input", "output","name-selected-region", "chr", "start", "end", "side-bp")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

# Load CpG counts
load(opt[["input"]])

# Filter the data for the given chromosome and region
CpG_count_region <- CpG_count[CpG_count$chr == opt[["chr"]] & CpG_count$start >= opt[["start"]] - opt[["side-bp"]] & CpG_count$start <= opt[["end"]] + opt[["side-bp"]], ]

# Save the CpG count filtered for the region of interest
save(CpG_count_region, file = opt[["output"]])
