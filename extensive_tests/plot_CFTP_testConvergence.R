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
  
  # Load summary stats nested replicates
  padded_n <- formatC(i, width = 2, format = "d", flag = "0")
  cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp_old.RData")
  load(cftp_file)
  cftp_sumstats_old <- summaryStats
  cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp.RData")
  load(cftp_file)
  cftp_sumstats_new <- summaryStats
  evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_old.RData")
  load(evol_file)
  evol_sumstats_old <- summaryStats
  evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, ".RData")
  load(evol_file)
  evol_sumstats_new <- summaryStats
  rm(summaryStats)
  
  # Load summary stats per replicate
  evol_sumstats_perRep <- list()
  cftp_sumstats_perRep <- list()
  
  for(rep in 1:10){
    
    file_name <- sprintf("summaryStats_CFTP_testConvergence%d_paramsID_%02d_rep_%02d.RData", test, i, rep)
    load(file_name)
    evol_sumstats_perRep[[rep]] <- summaryStats
    
    
    file_name <- sprintf("summaryStats_CFTP_testConvergence%d_paramsID_%02d_rep_%02d_cftp.RData", test, i, rep)
    load(file_name)
    cftp_sumstats_perRep[[rep]] <- summaryStats
  }
  
  ##### Frequency of P #####
  pdf(paste0("freqP_CFTP_testConvergence", test, "_paramsID_", padded_n, ".pdf"), width = 15, height = 10)
  
  # Set up layout for 2 rows, 4 columns, and space for a global title (oma parameter)
  par(mfrow = c(2, 4), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
  
  # Mean Island
  plot(evol_sumstats_old$meanFreqP_i, main = "Mean Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanFreqP_i, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanFreqP_i, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanFreqP_i, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanFreqP_i, main = "Mean Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanFreqP_i, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanFreqP_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanFreqP_i, col = colors[i], lty = 2, lwd = 2)
  }
  
  # Mean Non-Island
  plot(evol_sumstats_old$meanFreqP_ni, main = "Mean Non-Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanFreqP_ni, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanFreqP_ni, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanFreqP_ni, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanFreqP_ni, main = "Mean Non-Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanFreqP_ni, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanFreqP_ni, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanFreqP_ni, col = colors[i], lty = 2, lwd = 2)
  }
  
  # SD Island
  plot(evol_sumstats_old$sdFreqP_i, main = "SD Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$sdFreqP_i, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$sdFreqP_i, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$sdFreqP_i, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  # Plot summaryStats per replicate
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$sdFreqP_i, main = "SD Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$sdFreqP_i, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$sdFreqP_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$sdFreqP_i, col = colors[i], lty = 2, lwd = 2)
  }
  
  # SD Non-Island
  plot(evol_sumstats_old$sdFreqP_ni, main = "SD Non-Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$sdFreqP_ni, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$sdFreqP_ni, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$sdFreqP_ni, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  # Plot summaryStats per replicate
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$sdFreqP_ni, main = "SD Non-Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$sdFreqP_ni, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$sdFreqP_ni, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$sdFreqP_ni, col = colors[i], lty = 2, lwd = 2)
  }
  
  
  # Add global title for all plots after the first plot is created
  mtext("Frequency of partially-methylated sites", outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
  
  ##### Frequency of M #####
  pdf(paste0("freqM_CFTP_testConvergence", test, "_paramsID_", padded_n, ".pdf"), width = 15, height = 10)
  
  # Set up layout for 2 rows, 4 columns, and space for a global title (oma parameter)
  par(mfrow = c(2, 4), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
  
  # Mean Island
  plot(evol_sumstats_old$meanFreqM_i, main = "Mean Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanFreqM_i, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanFreqM_i, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanFreqM_i, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanFreqM_i, main = "Mean Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanFreqM_i, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanFreqM_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanFreqM_i, col = colors[i], lty = 2, lwd = 2)
  }
  
  # Mean Non-Island
  plot(evol_sumstats_old$meanFreqM_ni, main = "Mean Non-Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanFreqM_ni, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanFreqM_ni, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanFreqM_ni, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanFreqM_ni, main = "Mean Non-Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanFreqM_ni, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanFreqM_ni, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanFreqM_ni, col = colors[i], lty = 2, lwd = 2)
  }
  
  
  # SD Island
  plot(evol_sumstats_old$sdFreqM_i, main = "SD Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$sdFreqM_i, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$sdFreqM_i, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$sdFreqM_i, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  # Plot summaryStats per replicate
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$sdFreqM_i, main = "SD Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$sdFreqM_i, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$sdFreqM_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$sdFreqM_i, col = colors[i], lty = 2, lwd = 2)
  }
  
  # SD Non-Island
  plot(evol_sumstats_old$sdFreqM_ni, main = "SD Non-Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$sdFreqM_ni, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$sdFreqM_ni, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$sdFreqM_ni, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  # Plot summaryStats per replicate
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$sdFreqM_ni, main = "SD Non-Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$sdFreqM_ni, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$sdFreqM_ni, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$sdFreqM_ni, col = colors[i], lty = 2, lwd = 2)
  }
  
  
  # Add global title for all plots after the first plot is created
  mtext("Frequency of methylated sites", outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
  
  ##### Correlation #####
  
  pdf(paste0("Cor_Steepness_CFTP_testConvergence", test, "_paramsID_", padded_n, ".pdf"), width = 15, height = 10)
  
  # Set up layout for 2 rows, 4 columns, and space for a global title (oma parameter)
  par(mfrow = c(2, 4), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
  
  # Mean Island
  plot(evol_sumstats_old$meanCor_i, main = "Mean Cor. Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanCor_i, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanCor_i, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanCor_i, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanCor_i, main = "Mean Cor. Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanCor_i, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanCor_i, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanCor_i, col = colors[i], lty = 2, lwd = 2)
  }
  
  # Mean Non-Island
  plot(evol_sumstats_old$meanCor_ni, main = "Mean Cor. Non-Island (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanCor_ni, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanCor_ni, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanCor_ni, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanCor_ni, main = "Mean Cor. Non-Island (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanCor_ni, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanCor_ni, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanCor_ni, col = colors[i], lty = 2, lwd = 2)
  }
  
  
  
  ##### Steepness #####
  
  # Mean 
  plot(evol_sumstats_old$meanSteepness, main = "Mean Steepness (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$meanSteepness, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$meanSteepness, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$meanSteepness, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  transparent_colors <- sapply(colors, function(col) adjustcolor(col, alpha.f = 0.2))  # alpha.f is a transparency factor (0=fully transparent, 1=opaque)
  plot(evol_sumstats_perRep[[1]]$meanSteepness, main = "Mean Steepness (per rep.)", xlab = "Step Number", ylab = "Value", pch = 1, col = transparent_colors[1])
  abline(h = cftp_sumstats_perRep[[1]]$meanSteepness, col = colors[1], lty = 2, lwd = 2)
  for (i in 2:10){
    points(evol_sumstats_perRep[[i]]$meanSteepness, pch = i, col = transparent_colors[i])
    abline(h = cftp_sumstats_perRep[[i]]$meanSteepness, col = colors[i], lty = 2, lwd = 2)
  }
  
  # SD
  plot(evol_sumstats_old$sdSteepness, main = "SD Steepness (rep. average)", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
  points(evol_sumstats_new$sdSteepness, pch = 3, col = "blue")
  abline(h = cftp_sumstats_old$sdSteepness, col = "red", lty = 2, lwd = 2)
  abline(h = cftp_sumstats_new$sdSteepness, col = "blue", lty = 4, lwd = 2)
  legend("topright", 
         legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
         col = c("black", "red", "red", "black", "blue", "blue"), 
         lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
         pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
         ncol = 2,  # 2-column layout
         pt.cex = 1.5)  # Adjust point size if necessary
  
  
  # Add global title for all plots after the first plot is created
  mtext("Local correlations", outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
  
}

 ### OLD ###
if (FALSE){
  for (i in cases) {
    # Load summary stats nested replicates
    padded_n <- formatC(i, width = 2, format = "d", flag = "0")
    cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp_old.RData")
    load(cftp_file)
    cftp_sumstats_old <- summaryStats
    cftp_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_cftp.RData")
    load(cftp_file)
    cftp_sumstats_new <- summaryStats
    evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, "_old.RData")
    load(evol_file)
    evol_sumstats_old <- summaryStats
    evol_file <- paste0("summaryStats_nestedReps_CFTP_testConvergence", test, "_paramsID_", padded_n, ".RData")
    load(evol_file)
    evol_sumstats_new <- summaryStats
    rm(summaryStats)
    
    # Load summary stats per replicate
    evol_sumstats_perRep <- list()
    cftp_sumstats_perRep <- list()
    
    for(rep in 1:10){
      
      file_name <- sprintf("summaryStats_CFTP_testConvergence%d_paramsID_%02d_rep_%02d.RData", test, i, rep)
      load(file_name)
      evol_sumstats_perRep[[rep]] <- summaryStats
      
      
      file_name <- sprintf("summaryStats_CFTP_testConvergence%d_paramsID_%02d_rep_%02d_cftp.RData", test, i, rep)
      load(file_name)
      cftp_sumstats_perRep[[rep]] <- summaryStats
    }
    
    pdf(paste0("plot_CFTP_testConvergence", test, "_paramsID_", padded_n, ".pdf"), width = 15, height = 10)
    
    # Set up layout for 3 rows, 4 columns, and space for a global title (oma parameter)
    par(mfrow = c(3, 4), oma = c(0, 0, 3, 0))  # oma reserves space at the top for the title
    
    par(mfrow = c(1,2))
    # Plot each summary statistic, add axes labels and a legend
    plot(evol_sumstats_old$meanFreqP_i, main = "meanFreqP_i", xlab = "Step Number", ylab = "Value", pch = 1, col = "red")
    points(evol_sumstats_new$meanFreqP_i, pch = 3, col = "blue")
    abline(h = cftp_sumstats_old$meanFreqP_i, col = "red", lty = 2, lwd = 2)
    abline(h = cftp_sumstats_new$meanFreqP_i, col = "blue", lty = 4, lwd = 2)
    legend("topright", 
           legend = c("Old", "CFTP", "Evol", "New", "CFTP", "Evol"), 
           col = c("black", "red", "red", "black", "blue", "blue"), 
           lty = c(NA, 2, NA, NA, 4, NA),  # Line types for CFTP
           pch = c(NA, NA, 1, NA, NA, 3),  # Point types for Branch Evol
           ncol = 2,  # 2-column layout
           pt.cex = 1.5)  # Adjust point size if necessary
    
    # Plot summaryStats per replicate
    colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "yellow")
    plot(evol_sumstats_perRep[[1]]$meanFreqP_i, main = "meanFreqP_i", xlab = "Step Number", ylab = "Value", pch = 1, col = colors[1])
    for (i in 2:10){
      points(evol_sumstats_perRep[[i]]$meanFreqP_i, pch = i, col = colors[i])
      abline(h = cftp_sumstats_perRep[[i]]$meanFreqP_i, col = colors[i], lty = 2, lwd = 2)
    }
    par(mfrow = c(1,1))
    
    
    
    plot(evol_sumstats_old$meanFreqP_ni, main = "meanFreqP_ni", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanFreqP_ni, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$meanFreqM_i, main = "meanFreqM_i", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanFreqM_i, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$meanFreqM_ni, main = "meanFreqM_ni", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanFreqM_ni, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$sdFreqP_i, main = "sdFreqP_i", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$sdFreqP_i, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$sdFreqP_ni, main = "sdFreqP_ni", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$sdFreqP_ni, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$sdFreqM_i, main = "sdFreqM_i", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$sdFreqM_i, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$sdFreqM_ni, main = "sdFreqM_ni", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$sdFreqM_ni, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$meanCor_i, main = "meanCor_i", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanCor_i, col = "blue", lwd = 2)
    legend("bottomright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$meanCor_ni, main = "meanCor_ni", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanCor_ni, col = "blue", lwd = 2)
    legend("bottomright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$meanSteepness, main = "meanSteepness", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$meanSteepness, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    plot(evol_sumstats_old$sdSteepness, main = "sdSteepness", xlab = "Step Number", ylab = "Value", pch = 19)
    abline(h = cftp_sumstats_old$sdSteepness, col = "blue", lwd = 2)
    legend("topright", legend = c("CFTP", "Branch Evol"), col = c("blue", "black"), lty = c(1, NA), pch = c(NA, 19))
    
    # Add global title for all plots after the first plot is created
    mtext(paste("Parameter Combination:", i, ". Step length:", s), outer = TRUE, line = 1, cex = 1.5)
    
    dev.off()
    
  }
}


