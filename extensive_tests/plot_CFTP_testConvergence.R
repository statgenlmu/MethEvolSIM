test <- 2
cases <- c(2, 3, 4, 5, 7, 9, 10)
s <- 1 # step length used in the test

for (i in cases) {
  padded_n <- formatC(i, width = 2, format = "d", flag = "0")
  cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp.RData")
  load(cftp_file)
  cftp_sumstats <- summaryStats
  evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, ".RData")
  load(evol_file)
  evol_sumstats <- summaryStats
  rm(summaryStats)
  
  pdf(paste0("plot_CFTP_testConvergence", test, "_paramsID_", padded_n, ".pdf"), width = 15, height = 10)
  
  # Set up layout for 3 rows, 4 columns, and space for a global title (oma parameter)
  par(mfrow = c(3, 4), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
  
  # Individual plots with specific summary statistic titles
  plot(evol_sumstats$meanFreqP_i, main = "meanFreqP_i")
  abline(h = cftp_sumstats$meanFreqP_i)
  
  plot(evol_sumstats$meanFreqP_ni, main = "meanFreqP_ni")
  abline(h = cftp_sumstats$meanFreqP_ni)
  
  plot(evol_sumstats$meanFreqM_i, main = "meanFreqM_i")
  abline(h = cftp_sumstats$meanFreqM_i)
  
  plot(evol_sumstats$meanFreqM_ni, main = "meanFreqM_ni")
  abline(h = cftp_sumstats$meanFreqM_ni)
  
  plot(evol_sumstats$sdFreqP_i, main = "sdFreqP_i")
  abline(h = cftp_sumstats$sdFreqP_i)
  
  plot(evol_sumstats$sdFreqP_ni, main = "sdFreqP_ni")
  abline(h = cftp_sumstats$sdFreqP_ni)
  
  plot(evol_sumstats$sdFreqM_i, main = "sdFreqM_i")
  abline(h = cftp_sumstats$sdFreqM_i)
  
  plot(evol_sumstats$sdFreqM_ni, main = "sdFreqM_ni")
  abline(h = cftp_sumstats$sdFreqM_ni)
  
  plot(evol_sumstats$meanCor_i, main = "meanCor_i")
  abline(h = cftp_sumstats$meanCor_i)
  
  plot(evol_sumstats$meanCor_ni, main = "meanCor_ni")
  abline(h = cftp_sumstats$meanCor_ni)
  
  plot(evol_sumstats$meanSteepness, main = "meanSteepness")
  abline(h = cftp_sumstats$meanSteepness)
  
  plot(evol_sumstats$sdSteepness, main = "sdSteepness")
  abline(h = cftp_sumstats$sdSteepness)
  
  # Add global title for all plots after the first plot is created
  mtext(paste("Parameter Combination:", i, "Step length:", s), outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
}


