# Load necessary libraries
library(optparse)
source("functions_summaryStats.R") # Load functions to compute summary statistics

# Define the command-line options
option_list <- list(
  make_option(c("-d", "--data-dir"), type = "character", default = NULL,
              help = "Full path to the directory containing the .RData files", metavar = "character"),
  make_option(c("-f", "--design-file"), type = "character", default = NULL,
              help = "Full path to the simulation design file", metavar = "character"),
  make_option(c("-s", "--stats"), type = "character", default = "all",
              help = "Comma-separated list of summary statistics to compute (default: all). Options: meanFreqP_i, meanFreqP_ni, sdFreqP_i, sdFreqP_ni, meanFreqM_i, meanFreqM_ni, sdFreqM_i, sdFreqM_ni, meanFracMoverMU_i, meanFracMoverMU_ni, sdFracMoverMU_i, sdFracMoverMU_ni, FChangeCherry_i, FChangeCherry_ni, Fitch, Steepness, meanCor", metavar = "character"),
  make_option(c("-n", "--sample-n"), type = "integer", default = NULL,
              help = "Number of samples per file", metavar = "integer"),
  make_option(c("-p", "--pattern"), type = "character", default = NULL,
              help = ".RData files start name pattern", metavar = "character"),
  make_option(c("-u", "--update-file"), type = "character", default = NULL,
              help = "Full path to existing summary statistics file", metavar = "character")
)


# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt[["data-dir"]]) || is.null(opt[["design-file"]]) || is.null(opt[["sample-n"]]) || is.null(opt[["pattern"]])) {
  print_help(opt_parser)
  stop("All following arguments need to be provided: data-dir, design-file, sample-n, pattern")
}

# Convert the stats argument to a list
stats_to_compute <- strsplit(opt$stats, ",")[[1]]
stats_to_compute <- trimws(stats_to_compute)  # Remove any extra spaces
if (!is.null(opt[["update-file"]])){
  print("Updating summary statistics:")
  print(stats_to_compute)
} else {
  print("Computing summary statistics:")
  print(stats_to_compute)
}




##### SET INFORMATION STUDY DESIGN #############################################

# Load info simulation design: sampled_params, spatial_str, tree
load(opt[["design-file"]])

# Structural indeces of simulated CpG islands and non-islands
index_islands <- which(spatial_str$globalState=="U")
index_nonislands <- which(spatial_str$globalState=="M")

# Number of CpGs
if(length(unique(spatial_str$n)) == 1){
  n_CpG <- unique(spatial_str$n)
}else{
  print("Different number of CpGs per simulated structure")
}

##### IMPORT DATA AND EXTRACT SUMMARY STATISTICS ###############################

# List simulation output files
RData_files <- list.files(opt[["data-dir"]], pattern = paste0("^", opt[["pattern"]], "_cftp\\.RData$"), full.names = TRUE)



# Set the number of simulations to process
n_sim <- length(RData_files)

if (!is.null(opt[["update-file"]])){
  # Load existing summary statistics file to update
  load(opt[["update-file"]])
  
} else {
  # Initialize the vectors to store the summary statistics
  meanFreqP_i <- rep(NA, n_sim)
  meanFreqP_ni <- rep(NA, n_sim)
  sdFreqP_i <- rep(NA, n_sim)
  sdFreqP_ni <- rep(NA, n_sim)
  meanFreqM_i <- rep(NA, n_sim)
  meanFreqM_ni <- rep(NA, n_sim)
  sdFreqM_i <- rep(NA, n_sim)
  sdFreqM_ni <- rep(NA, n_sim)
  meanFracMoverMU_i <- rep(NA, n_sim)
  meanFracMoverMU_ni <- rep(NA, n_sim)
  sdFracMoverMU_i <- rep(NA, n_sim)
  sdFracMoverMU_ni <- rep(NA, n_sim)
  
  meanFChangeCherry_rate_i <- rep(NA, n_sim)
  meanFChangeCherry_scale_i <- rep(NA, n_sim)
  meanFChangeCherry_rate_ni <- rep(NA, n_sim)
  meanFChangeCherry_scale_ni <- rep(NA, n_sim)
  
  meanFitch <- rep(NA, n_sim)
  sdFitch <- rep(NA, n_sim)
  
  meanSteepness <- rep(NA, n_sim)
  sdSteepness <- rep(NA, n_sim)
  
  meanCor_i <- rep(NA, n_sim)
  meanCor_ni <- rep(NA, n_sim)
  #meanCov_i <- rep(NA, n_sim)
  #meanCov_ni <- rep(NA, n_sim)  
} 




error_log <- list()  # To store error messages
print(paste("Computing summary statistics for", n_sim, "simulations."))
print("Files to process:")
print(RData_files)
for(sim in 1:n_sim){
  print(paste("processing simulation:", sim))
  load(RData_files[sim])
  
  # Compute summary Statistics for the frequency of P
  if ( "all" %in% stats_to_compute || "meanFreqP_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFreqP_i"] <- get_islandMeanFreqP(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFreqP_i[sim] <- get_islandMeanFreqP(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "meanFreqP_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFreqP_ni"] <- get_nonislandMeanFreqP(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFreqP_ni[sim] <- get_nonislandMeanFreqP(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFreqP_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFreqP_i"] <- get_islandSDFreqP(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFreqP_i[sim] <- get_islandSDFreqP(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFreqP_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFreqP_ni"] <- get_nonislandSDFreqP(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFreqP_ni[sim] <- get_nonislandSDFreqP(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  
  # Compute summary statistics for the frequency of M
  if ( "all" %in% stats_to_compute || "meanFreqM_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFreqM_i"] <- get_islandMeanFreqM(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFreqM_i[sim] <- get_islandMeanFreqM(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "meanFreqM_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFreqM_ni"] <- get_nonislandMeanFreqM(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFreqM_ni[sim] <- get_nonislandMeanFreqM(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFreqM_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFreqM_i"] <- get_islandSDFreqM(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFreqM_i[sim] <- get_islandSDFreqM(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFreqM_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFreqM_ni"] <- get_nonislandSDFreqM(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFreqM_ni[sim] <- get_nonislandSDFreqM(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  
  # Compute summary statistics for the fraction of M over MU
  if ( "all" %in% stats_to_compute || "meanFracMoverMU_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFracMoverMU_i"] <- get_islandMeanFracMoverMU(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFracMoverMU_i[sim] <- get_islandMeanFracMoverMU(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "meanFracMoverMU_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFracMoverMU_ni"] <- get_nonislandMeanFracMoverMU(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanFracMoverMU_ni[sim] <- get_nonislandMeanFracMoverMU(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFracMoverMU_i" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFracMoverMU_i"] <- get_islandSDFracMoverMU(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFracMoverMU_i[sim] <- get_islandSDFracMoverMU(index_islands = index_islands, data = data, sample_n = opt[["sample-n"]])
    }
  }
  if ( "all" %in% stats_to_compute || "sdFracMoverMU_ni" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"sdFracMoverMU_ni"] <- get_nonislandSDFracMoverMU(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    } else {
      sdFracMoverMU_ni[sim] <- get_nonislandSDFracMoverMU(index_nonislands = index_nonislands, data = data, sample_n = opt[["sample-n"]])
    }
  }

  
  # Compute summary statistics of tree cherry comparisons
  if ( "all" %in% stats_to_compute || "FChangeCherry_i" %in% stats_to_compute) {
    FChange_cherryData <- get_FChange_cherryData(tree = tree, data = data, n_CpG = n_CpG)
    
    meanFChangeCherry_fit_i <- tryCatch({
      fit_meanFChange_cherry_i(FChange_cherryData = FChange_cherryData, index_islands = index_islands)
    }, error = function(e) {
      error_msg <- paste("Simulation", sim, "- Error in fit_meanFChange_cherry_i:", e$message)
      error_log[[length(error_log) + 1]] <- error_msg
      list(rate = NA, scale = NA)
    })
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFChangeCherry_rate_i"] <- meanFChangeCherry_fit_i$rate
      summaryStats[sim,"meanFChangeCherry_scale_i"] <- meanFChangeCherry_fit_i$scale
    } else {
      meanFChangeCherry_rate_i[sim] <- meanFChangeCherry_fit_i$rate
      meanFChangeCherry_scale_i[sim] <- meanFChangeCherry_fit_i$scale
    }
  }
  if ( "all" %in% stats_to_compute || "FChangeCherry_ni" %in% stats_to_compute) {
    meanFChangeCherry_fit_ni <- tryCatch({
      fit_meanFChange_cherry_ni(FChange_cherryData = FChange_cherryData, index_nonislands = index_nonislands)
    }, error = function(e) {
      error_msg <- paste("Simulation", sim, "- Error in fit_meanFChange_cherry_ni:", e$message)
      error_log[[length(error_log) + 1]] <- error_msg
      list(rate = NA, scale = NA)
    })
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFChangeCherry_rate_ni"] <- meanFChangeCherry_fit_ni$rate
      summaryStats[sim,"meanFChangeCherry_scale_ni"] <- meanFChangeCherry_fit_ni$scale
    } else {
      meanFChangeCherry_rate_ni[sim] <- meanFChangeCherry_fit_ni$rate
      meanFChangeCherry_scale_ni[sim] <- meanFChangeCherry_fit_ni$scale
    }
  }
  

  # Compute summary statistics with Fitch algorithm for global state of islands
  if ( "all" %in% stats_to_compute || "Fitch" %in% stats_to_compute) {
    fitchData <- computeFitch_RegionGlbSt(index_islands = index_islands, data = data, tree = tree, u_threshold = 0.2, m_threshold = 0.8)
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanFitch"] <- mean(fitchData)
      summaryStats[sim,"sdFitch"] <- sd(fitchData)
    } else {
      meanFitch[sim] <- mean(fitchData)
      sdFitch[sim] <- sd(fitchData)
    }
  }
  
  ## Compute summary statistics for the transitions of methylation frequency between differentially methylated regions
  if ( "all" %in% stats_to_compute || "Steepness" %in% stats_to_compute) {
    transition_params <- tryCatch({
      fit_MethTrans(data = data, threshold = 0.5, minRepresentation = 20, subset_CpG_n = 30)
    }, error = function(e) {
      error_msg <- paste("Simulation", sim, "- Error in fit_MethTrans:", e$message)
      error_log[[length(error_log) + 1]] <- error_msg
      list(meanSteepness = NA, sdSteepness = NA)
    })
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanSteepness"] <- transition_params$meanSteepness
      summaryStats[sim,"sdSteepness"] <- transition_params$sdSteepness
    } else {
      meanSteepness[sim] <- transition_params$meanSteepness
      sdSteepness[sim] <- transition_params$sdSteepness
    }
  }
  
  
  ## Correlation/Covariance within structure type
  if ( "all" %in% stats_to_compute || "meanCor" %in% stats_to_compute) {
    if (!is.null(opt[["update-file"]])){
      summaryStats[sim,"meanCor_i"] <- compute_meanCor_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
      summaryStats[sim,"meanCor_ni"] <- compute_meanCor_i(index_nonislands = index_nonislands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
    } else {
      meanCor_i[sim] <- compute_meanCor_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
      meanCor_ni[sim] <- compute_meanCor_ni(index_nonislands = index_nonislands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
    }
  }
  #if ( "all" %in% stats_to_compute || "meanCov" %in% stats_to_compute) {
  #  if (!is.null(opt[["update-file"]])){
    #summaryStats[sim,"meanCov_i"] <- <- compute_meanCov_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
    #ummaryStats[sim,"meanCov_i"] <- <- compute_meanCov_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
  #} else {
    #meanCov_i[sim] <- compute_meanCov_i(index_islands = index_islands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
    #meanCov_ni[sim] <- compute_meanCov_ni(index_nonislands = index_nonislands, minN_CpG = 80, shore_length = 10, data = data, sample_n = opt[["sample-n"]])
  #}
  #}
}  

# Save the error log if there are any errors
if (length(error_log) > 0) {
  print("Errors during summary statistics computation printed in file error_log.txt")
  writeLines(unlist(error_log), "error_log.txt")
}  
  
if (!is.null(opt[["update-file"]])){
  # Save the updated dataframe
  print("Finished updating summary statitistics")
} else {
  # Save the vectors as dataframe
  summaryStats <- data.frame(meanFreqP_i = meanFreqP_i,
                             meanFreqP_ni = meanFreqP_ni,
                             sdFreqP_i = sdFreqP_i,
                             sdFreqP_ni = sdFreqP_ni,
                             meanFreqM_i = meanFreqM_i,
                             meanFreqM_ni = meanFreqM_ni,
                             sdFreqM_i = sdFreqM_i,
                             sdFreqM_ni = sdFreqM_ni,
                             meanFracMoverMU_i = meanFracMoverMU_i,
                             meanFracMoverMU_ni = meanFracMoverMU_ni,
                             sdFracMoverMU_i = sdFracMoverMU_i,
                             sdFracMoverMU_ni = sdFracMoverMU_ni,
                             meanFChangeCherry_rate_i = meanFChangeCherry_rate_i,
                             meanFChangeCherry_scale_i = meanFChangeCherry_scale_i,
                             meanFChangeCherry_rate_ni = meanFChangeCherry_rate_ni,
                             meanFChangeCherry_scale_ni = meanFChangeCherry_scale_ni,
                             meanFitch = meanFitch,
                             sdFitch = sdFitch,
                             meanSteepness = meanSteepness,
                             sdSteepness = sdSteepness,
                             meanCor_i = meanCor_i,
                             meanCor_ni = meanCor_ni)
  print("Finished generating summary statistics")
} 

# Set output name and save
out_name <- paste0("summaryStats_",opt[["pattern"]], "_cftp.RData")
save(summaryStats, file = out_name)
print(paste("Finished processing. Generated file:", out_name, "under directory:", getwd()))  
  
  
  
  
  




