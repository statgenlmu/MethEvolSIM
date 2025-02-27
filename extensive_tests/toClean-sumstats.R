

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
