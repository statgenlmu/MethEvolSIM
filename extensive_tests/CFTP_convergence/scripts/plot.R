library(optparse)


# Define the command-line options
option_list <- list(
  make_option("--dir", type = "character", default = NULL,
              help = "Full path to save the output file", metavar = "character"),
  make_option(c("-f", "--design-file-name"), type = "character", default = NULL,
              help = "Name of the simulation design file", metavar = "character"),
  make_option("--n-sim", type = "integer", default = NULL,
              help = "Number of simulations", metavar = "integer"),
  make_option(c("-r", "--replicate-n"), type = "integer", default = NULL,
              help = "number of data replicates to simulate", metavar = "integer")
)

# Parse the arguments

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Get the names of required options (you can update this based on what's mandatory)
required_options <- c("dir", "design-file-name", "n-sim", "replicate-n")

# Check that all required options are not NULL
missing_options <- required_options[sapply(required_options, function(x) is.null(opt[[x]]))]

if (length(missing_options) > 0) {
  print_help(opt_parser)
  stop(paste("The following arguments need to be provided:", paste(missing_options, collapse = ", ")))
}

# Assign the options to variables
dir <- opt[["dir"]]
design_file <- opt[["design-file-name"]]
n_sim <- opt[["n-sim"]]
n_rep <- opt[["replicate-n"]]

# Load the design file
load(file.path(dir, design_file))

#sampled_params$iota[order(sampled_params$iota, decreasing = TRUE)]

# Set a list to save the summary statistics for the replicates of each parameter combination
branchEvol_sumstats <- list()

# Set the pad number for parameter combination and replicate number
params_pad_n <- nchar(as.character(n_sim))
replicate_pad_n = nchar(as.character(n_rep))

# Set the colors to use in the plot
colors <- rainbow(n_rep)
transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)


for(index_params in 1:n_sim){
  padded_index_params <- formatC(index_params, width = params_pad_n, format = "d", flag = "0")
  
  for(rep in 1:n_rep){
    padded_replicate_n <- formatC(rep, width = replicate_pad_n, format = "d", flag = "0")
    
    # Load the summary statistics of branch_evol and save them
    load(file.path(dir, paste0("summaryStats_CFTP_testConvergence_paramsID_", padded_index_params, "_rep_", padded_replicate_n, ".RData")))
    branchEvol_sumstats[[rep]] <- summaryStats
  }
  # Load the ksummary statistics of the cftp method and save them
  load(file.path(dir, paste0("summaryStats_CFTP_testConvergence_paramsID_", padded_index_params, "_cftp.RData")))
  cftp_sumstats <- summaryStats
  rm(summaryStats)
  
  pdf(paste0("Figures/CFTP_testConvergence_paramsID_", padded_index_params,".pdf"), width = 8, height = 4)
  
  # Set up layout for 1 rows, 2 columns, and space for a global title (oma parameter)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
  
  # Plot island correlations
  plot(branchEvol_sumstats[[1]]$meanCor_i, main = "Mean Cor. Island", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats$meanCor_i[1], col = colors[1], lty = 2, lwd = 2)
  for (i in 2:n_rep){
    points(branchEvol_sumstats[[i]]$meanCor_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats$meanCor_i[i], col = colors[i], lty = 2, lwd = 2)
  }
  
  # Plot non-island correlations
  plot(branchEvol_sumstats[[1]]$meanCor_ni, main = "Mean Cor. Non-Island", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats$meanCor_i[1], col = colors[1], lty = 2, lwd = 2)
  for (i in 2:n_rep){
    points(branchEvol_sumstats[[i]]$meanCor_in, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats$meanCor_i[i], col = colors[i], lty = 2, lwd = 2)
  }
  
  # Add global title for all plots after the first plot is created
  mtext("Local correlations", outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
}


