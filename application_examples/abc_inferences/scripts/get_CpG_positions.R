library(optparse)
library(Biostrings)

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input CpG_count.RData path", metavar = "character"),
  make_option(c("--reference-genome"), type = "character", default = NULL,
              help = "Relative path to the reference genome file (.fa)", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output relative path with file name", metavar = "character")
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


# Load Reference Genome and CpG_count
genome <- readDNAStringSet(opt[["reference-genome"]], format = "fasta")
load(opt[["input"]])

# Define a function to get CpG positions in a given sequence
get_CpG_positions <- function(seq, start_pos) {
  cpg_sites <- start(Biostrings::matchPattern("CG", seq))  # Find "CG" occurrences
  return(cpg_sites + start_pos - 1)  # Convert to genome-wide positions
}

# Extract only the chromosome part from reference genome names
genome_chr_names <- sapply(names(genome), function(x) strsplit(x, " ")[[1]][1])

# Filter only the regions with at least one CpG site
CpG_count <- CpG_count[CpG_count$CpG_count > 0,]

# Initialize a list to store results
cpg_positions_list <- list()

# Initialize a variable to keep track of current chromosome
current_chr <- "no_chromosome"

# Loop over each region in CpG_count
for (i in 1:nrow(CpG_count)) {
  region <- CpG_count[i, ]
  chr_name <- region$chr
  
  # Extract chromosome sequence
  if (chr_name != current_chr){
    
    current_chr <- chr_name
    print(paste("Processing chromosome", current_chr))
    
    # Check if CpG_count chromosome names match genome_chr_names
    chr_index <- which(genome_chr_names == current_chr)
    
    if (length(chr_index) > 0) {
      seq_chr <- genome[[chr_index]]  # Extract full chromosome sequence
    } else {
      print(paste("Chromosome", current_chr, "not found in genome!"))
    }
  }
  
  # Extract region sequence
  seq_region <- subseq(seq_chr, start = region$start, end = region$end)  # Extract region
  
  # Get CpG positions
  cpg_positions <- get_CpG_positions(seq_region, region$start)
  
  # Store results
  cpg_positions_list[[i]] <- data.frame(
    chr = chr_name,
    region_type = region$str,
    CpG_position = cpg_positions
  )
}

# Combine results into a single data frame and save into output file
cpg_positions <- do.call(rbind, cpg_positions_list)
save(cpg_positions, file = opt[["output"]])

