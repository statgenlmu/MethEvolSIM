library(optparse)

option_list <- list(
  make_option(c("--filtered-CpGcount-region"), type = "character", default = NULL,
              help = "CpG_count_region.RData relative path", metavar = "character"),
  make_option(c("--CpG-positions"), type = "character", default = NULL,
              help = "CpG_positions.RData relative path", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output relative path with file name", metavar = "character")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("filtered-CpGcount-region", "CpG-positions", "output")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}


# Load the CpG count for the genomic region of interest and the CpG positions 
# relative to the reference genome (cytosine, relative to the positive strand)
load(opt[["filtered-CpGcount-region"]])
load(opt[["CpG-positions"]])

# Initialize a list to store results
filtered_positions <- list(length = nrow(CpG_count_region))

# Loop over each island and non-island region in the genomic region of interest
for (i in 1:nrow(CpG_count_region)) {
  
  # Extract the relevant information for the current island / non-island region
  region <- CpG_count_region[i,]
  chr <- region$chr
  start <- region$start
  end <- region$end
  
  # Extract the CpG positions within the current region
  filtered_positions[[i]] <- cpg_positions[cpg_positions$chr == chr & cpg_positions$CpG_position >= start & cpg_positions$CpG_position <= end,]
  
  # Add the positions for region start and end
  filtered_positions[[i]]$region_start <- region$start
  filtered_positions[[i]]$region_end <- region$end
  
  # Check the number of filtered positions corresponds to the CpG count for that region
  if(nrow(filtered_positions[[i]]) != region$CpG_count){
    print("Numper of filtered positions differs from region CpG_count in region:")
    print(region)
  }
  
  # Add the CpG count for that region
  filtered_positions[[i]]$region_CpG_count <- region$CpG_count
  
  # Add the CpG relative index within the region
  filtered_positions[[i]] <- filtered_positions[[i]][order(filtered_positions[[i]]$CpG_position), ]
  filtered_positions[[i]]$relative_index <- seq_len(nrow(filtered_positions[[i]]))
}

# Combine results into a single data frame and save into output file
filtered_positions <- do.call(rbind, filtered_positions)
save(filtered_positions, file = opt[["output"]])


