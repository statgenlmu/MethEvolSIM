test <- 2
if (test == 1){
  s <- 0.1 # step length used in the test
  n <- 500 # number of steps used in the test
} else if (test == 2){
  s <- 1 # step length used in the test
  n <- 250 # number of steps used in the test
}
cases <- c(2, 3, 4, 5, 7, 9, 10)


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
  
  # Plot each summary statistic, add axes labels and a legend
  plot(evol_sumstats$meanFreqP_i, main = "meanFreqP_i", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanFreqP_i, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanFreqP_ni, main = "meanFreqP_ni", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanFreqP_ni, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanFreqM_i, main = "meanFreqM_i", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanFreqM_i, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanFreqM_ni, main = "meanFreqM_ni", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanFreqM_ni, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$sdFreqP_i, main = "sdFreqP_i", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$sdFreqP_i, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$sdFreqP_ni, main = "sdFreqP_ni", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$sdFreqP_ni, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$sdFreqM_i, main = "sdFreqM_i", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$sdFreqM_i, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$sdFreqM_ni, main = "sdFreqM_ni", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$sdFreqM_ni, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanCor_i, main = "meanCor_i", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanCor_i, col = "blue", lwd = 2)
  legend("bottomright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanCor_ni, main = "meanCor_ni", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanCor_ni, col = "blue", lwd = 2)
  legend("bottomright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$meanSteepness, main = "meanSteepness", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$meanSteepness, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  plot(evol_sumstats$sdSteepness, main = "sdSteepness", xlab = "Step Number", ylab = "Value", pch = 19)
  abline(h = cftp_sumstats$sdSteepness, col = "blue", lwd = 2)
  legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
  
  # Add global title for all plots after the first plot is created
  mtext(paste("Parameter Combination:", i, ". Step length:", s), outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
  
}

## check convergence meanSteepness
for (i in cases){
  print(paste("Test:", test, ". Params comb:", i, ". Step length:", s, ". Step number:", n))
  padded_n <- formatC(i, width = 2, format = "d", flag = "0")
  cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp.RData")
  load(cftp_file)
  cftp_sumstats <- summaryStats
  evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, ".RData")
  load(evol_file)
  evol_sumstats <- summaryStats
  rm(summaryStats)
  print(paste("cftp_sumstats$meanSteepness:", cftp_sumstats$meanSteepness))
  if(n == 250){
    print(paste("mean(evol_sumstats$meanSteepness[150:250])", mean(evol_sumstats$meanSteepness[150:250])))
    print(paste("sd(evol_sumstats$meanSteepness[150:250])", sd(evol_sumstats$meanSteepness[150:250])))
  } else if (n == 500){
    print(paste("mean(evol_sumstats$meanSteepness[150:500])", mean(evol_sumstats$meanSteepness[150:500])))
    print(paste("sd(evol_sumstats$meanSteepness[150:500])", sd(evol_sumstats$meanSteepness[150:500])))
  }
  
}

## Possibilities for not complete convergence:
# 1. cftp$meanSteepness is averaged over 10 replicates, evol too but at many time steps
# 2. 
load("/home/sara/methylation/phd_project_saracv/simulation_studies/abc_designSIM.RData")
sampled_params[7,] # but with mu = 0



