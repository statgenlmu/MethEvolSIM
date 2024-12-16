library(optparse)
library(Biostrings) # For reading the genome sequence

option_list <- list(
  make_option(c("--CGI-annotation"), type = "character", default = NULL,
              help = "Path to the CGI annotation design file (.txt)", metavar = "character"),
  make_option(c("--reference-genome"), type = "character", default = NULL,
              help = "Path to the reference genome file (.fa)", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output path", metavar = "character")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("CGI-annotation", "reference-genome", "output")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]


# Step 1: Load Reference Genome and CGI Annotation
cgi_data <- read.table(opt[["CGI-annotation"]], header = TRUE, sep = "\t")
genome <- readDNAStringSet(opt[["reference-genome"]], format = "fasta")

# Step 2: Define a function to compute CpG counts in a region
count_CpG <- function(sequence) {
  Biostrings::countPattern("CG", sequence)
}

# Step 3: Generate a Count Data Frame

# Initialize the dataframe
CpG_count <- data.frame(
  chr = character(),
  str = factor(),
  start = numeric(),
  end = numeric(),
  CpG_count = numeric(),
  metadata = character(),
  stringsAsFactors = FALSE
)

# Iterate over the chromosomes in the genome
for (chrom in names(genome)) {
  # Get chromosome sequence
  chr_seq <- genome[[chrom]]
  chr_len <- length(chr_seq)
  
  # Filter CGI data for this chromosome
  cleaned_chrom <- strsplit(chrom, " ")[[1]][1]
  cgi_chr <- cgi_data[cgi_data$chr == cleaned_chrom, ]
  
  # Check if there are any CGIs for this chromosome
  if (nrow(cgi_chr) == 0) {
    next  # Skip to the next chromosome if no data is available for this one
  }
  
  # Extract metadata
  metadata <- paste(strsplit(chrom, " ")[[1]][-1], collapse = " ")
  
  # Initialize previous end position (to calculate non-islands)
  prev_end <- 0
  
  # Process each CGI region
  for (i in 1:nrow(cgi_chr)) {
    island_start <- cgi_chr$start[i]
    island_end <- cgi_chr$end[i]
    island_CpG_count <- cgi_chr$CpGcount[i]
    
    # Define the non-island before this island
    if (island_start > prev_end + 1) {
      non_island_start <- prev_end + 1
      non_island_end <- island_start - 1
      non_island_seq <- subseq(chr_seq, non_island_start, non_island_end)
      non_island_CpG_count <- count_CpG(non_island_seq)
      
      # Add non-island row to restructured output
      CpG_count <- rbind(CpG_count, data.frame(
        chr = cleaned_chrom,
        str = factor(0, levels = c(0, 1), labels = c("Non-island", "Island")),
        start = non_island_start,
        end = non_island_end,
        CpG_count = non_island_CpG_count,
        metadata = metadata
      ))
    }
    
    # Add island row to restructured output
    CpG_count <- rbind(CpG_count, data.frame(
      chr = cleaned_chrom,
      str = factor(1, levels = c(0, 1), labels = c("Non-island", "Island")),
      start = island_start,
      end = island_end,
      CpG_count = island_CpG_count,
      metadata = metadata
    ))
    
    # Update previous end
    prev_end <- island_end
  }
  
  # Handle non-island region after the last island
  if (prev_end < chr_len) {
    non_island_start <- prev_end + 1
    non_island_end <- chr_len
    non_island_seq <- subseq(chr_seq, non_island_start, non_island_end)
    non_island_CpG_count <- count_CpG(non_island_seq)
    
    # Add final non-island row to restructured output
    CpG_count <- rbind(CpG_count, data.frame(
      chr = cleaned_chrom,
      str = factor(0, levels = c(0, 1), labels = c("Non-island", "Island")),
      start = non_island_start,
      end = non_island_end,
      CpG_count = non_island_CpG_count,
      metadata = metadata
    ))
  }
}

# Sort the output by chromosome and start position
CpG_count <- CpG_count[order(CpG_count$chr, CpG_count$start), ]

# Save the output
save(CpG_count, file = opt[["output"]])
