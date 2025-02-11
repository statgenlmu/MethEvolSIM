
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






#undebug(countSites_cherryMethDiff)
#undebug(get_cherryDist)

##TODO: n_CpG not given but calculated. Each structure is assumed to have the same number of sites at all tips
get_siteFChange_cherry_perMethDiffType <- function(tree, data, sample_n){
  if(sample_n < 2){stop("Minimum number of tips/samples needed: 2")}
  cherries <- get_cherryDist(tree = tree, sample_n = sample_n)
  siteFChange_cherry_perMethDiffType <- countSites_cherryMethDiff(cherries = cherries, data = data)
  str_n <- length(data[[1]]) # set the number of structures
  sites_n <- numeric(length = str_n)
  for(str in 1:str_n) sites_n[str] <- length(data[[1]][[str]])
  # Each value is a count number of full methylation changes _f or changes in one strand _h.
  # Each structure has a number of CpGs. Then, if we wanna count the proportion
  # of sites, we divide each value by the number of CpGs. 
  # Explicitly match sites_n to the columns of siteFChange_cherry_perMethDiffType
  sites_n_repeated <- rep(sites_n, each = 2)
  siteFChange_cherry_perMethDiffType[,3:ncol(siteFChange_cherry_perMethDiffType)] <- sweep(siteFChange_cherry_perMethDiffType[,3:ncol(siteFChange_cherry_perMethDiffType)], 2, sites_n_repeated, "/")
  siteFChange_cherry_perMethDiffType
}

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
