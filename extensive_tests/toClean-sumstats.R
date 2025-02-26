
index_islands <- c(1,3)
index_nonislands <- c(2,4)
# Simulated data has all positions
# But empirical data may have differences in coverage
# E.g. if in island one positions 1,3,5,7 are CpGs
# Empirical data may have in first_tip 1,5,7
# And in second tip 1,3,5,7
# So that only sites 1,5,7 can be compared
# Previous step of filtering data using a list containing for each island and non-island region a vector of the
# relative indices of the CpGs that are to be compared.




get_siteFChange_cherry <- function(tree, data, sample_n){
  if(sample_n < 2){stop("Minimum number of tips/samples needed: 2")}
  siteFChange_cherry_perMethDiffType <- get_siteFChange_cherry_perMethDiffType(tree = tree, data = data, sample_n = sample_n)
  # Identify the unique structures in the column names
  structures <- unique(sub("_[fh]$", "", names(siteFChange_cherry_perMethDiffType)[-c(1, 2)]))
  # Create the new dataframe
  siteFChange_cherry <- data.frame(
    tips = siteFChange_cherry_perMethDiffType$tips,
    dist = siteFChange_cherry_perMethDiffType$dist
  )
  # Add the summed columns
  for (str in structures){
    siteFChange_cherry[[str]] <- rowSums(siteFChange_cherry_perMethDiffType[, grep(paste0("^", str, "_"), names(siteFChange_cherry_perMethDiffType))])
  }
  siteFChange_cherry
}

## For something like before:
# There are 2 columns per structure in siteFChange_cherry_perMethDiffType
# _f represents the frequency of full changes (both DNA strands)
# _h represents the frequency of half changes (one DNA strand)
#siteFChange_cherry_perMethDiffType[,grep("h", colnames(siteFChange_cherry_perMethDiffType))] <- siteFChange_cherry_perMethDiffType[,grep("h", colnames(siteFChange_cherry_perMethDiffType))]/2

#index islands: vector

SDandMean_SiteFChange_cherry <- function(siteFChange_cherry, index_islands, index_nonislands){
  siteFChange <- siteFChange_cherry[-c(1,2)]
  if(length(index_nonislands)>0){
    mean_NI <- rowMeans(as.data.frame(siteFChange[,index_nonislands]))
    sd_NI <- apply(as.data.frame(siteFChange[,index_nonislands]), 1, sd)
  } else {
    mean_NI = sd_NI <- NA
  }
  if(length(index_islands)>0){
    mean_I <- rowMeans(as.data.frame(siteFChange[,index_islands]))
    sd_I <- apply(as.data.frame(siteFChange[,index_islands]), 1, sd)
  } else {
    mean_I = sd_I <- NA
  }
  SDandMean_SiteFChange_cherry  <- data.frame(
    tips = siteFChange_cherry$tips,
    dist = siteFChange_cherry$dist,
    mean_I = mean_I,
    sd_I = sd_I,
    mean_NI = mean_NI,
    sd_NI = sd_NI
  )
  SDandMean_SiteFChange_cherry
}
#undebug(SDandMean_SiteFChange_cherry)
##TODO: Add index_nonislands
