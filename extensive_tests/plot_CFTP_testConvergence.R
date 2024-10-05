test <- 1
cases <- c(2,3,4,5,7,9,10)

for (i in cases){
  padded_n <- formatC(i, width = 2, format = "d", flag = "0")
  cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp.RData")
  load(cftp_file)
  cftp_sumstats <- summaryStats
  evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, ".RData")
  load(evol_file)
  evol_sumstats <- summaryStats
  rm(summaryStats)
  pdf(paste0("plot_CFTP_testConvergence", test, "_paramsID_", padded_n,".pdf"), width = 15, height = 6)
  par(mfrow = c(2,5))
  plot(evol_sumstats$meanFreqP_i, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanFreqP_i)
  plot(evol_sumstats$meanFreqP_ni, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanFreqP_ni)
  plot(evol_sumstats$meanFreqM_i, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanFreqM_i)
  plot(evol_sumstats$meanFreqM_ni, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanFreqM_ni)
  plot(evol_sumstats$sdFreqP_i, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$sdFreqP_i)
  plot(evol_sumstats$sdFreqP_ni, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$sdFreqP_ni)
  plot(evol_sumstats$sdFreqM_i, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$sdFreqM_i)
  plot(evol_sumstats$sdFreqM_ni, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$sdFreqM_ni)
  plot(evol_sumstats$meanCor_i, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanCor_i)
  plot(evol_sumstats$meanCor_ni, main = paste("Parameter Combination:", i))
  abline(h = cftp_sumstats$meanCor_ni)
  dev.off()
}  

