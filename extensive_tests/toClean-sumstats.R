

#### #### #### Tree cherries comparisons #### #### ####

find_cherries <- function(tree) {
  ## tree must be an ape tree whose tip labels are numbers.
  ## returns a list of vectors (a, b, x), where  (a, b) are a cherry of tree,
  ## that is a pair of leaves that are sisters in the tree,
  ## and x is their distance, that is, sum of the branch lengths between them 
  outlist <- list()
  d <- cophenetic.phylo(tree)
  n <- nrow(d)
  for(i in 1:n) {
    j <- which(d[i,]==min(d[i, -i]))
    if(length(j)==1) {
      k <- which(d[j,]==min(d[j, -j]))
      if(length(k)==1) {
        if(k==i & i<j) {
          ## outlist[[length(outlist)+1]] <- c(as.integer(rownames(d)[i]), as.integer(colnames(d)[j]), d[i,j])
          outlist[[length(outlist)+1]] <- c(i, j, d[i,j])
        }
      }
    }
  }
  outlist
}



## TODO: The tree needs to
## have the tips with numeric labels. ordered from left to as ((1:1,2:1):5,3:6); 
## For that check previous (old) function 
## tree needs to be an object of ape's phylo class or a character string in
## parenthetic format known as the Newick or New Hampshire format
get_cherryDist <- function(tree, sample_n, testing = FALSE){
  if(sample_n < 2){stop("Minimum number of tips/samples needed: 2")}
  if(class(tree) != "phylo") tree <- ape::read.tree(text = tree)
  # compute the pairwise distances between the tips from a phylogenetic tree
  dist <- ape::cophenetic.phylo(tree)
  # get the number of tips
  n_tips <- nrow(dist)
  # set list to store the cherry info
  cherry_list <- list()
  for (tip in 1:n_tips){
    # find the closest tip
    tip_a <- which(dist[tip,] == min(dist[tip, -tip]))
    # if there is one single closest tip
    if(length(tip_a)==1){
      # find the closes tip 
      tip_b <- which(dist[tip_a,] == min(dist[tip_a, -tip_a]))
      # if there is one single closest tip
      if (length(tip_b)==1){
        # if both are mutually closest tips and
        # tip < tip_a is used to avoid identifying the cherry twice (dist is symmetric)
        if(tip_b == tip & tip < tip_a){
          cherry_list[[length(cherry_list)+1]] <- c(tip, tip_a, dist[tip,tip_a])
        }
      }
    }
  }
  if(testing){
    list(tree = tree,
         cherry_list = cherry_list)
  } else {
    cherry_list
  }
  
}


#undebug(get_cherryDist)
## TODO: Assumes that each str has the same number of sites at the different tips
## Think about how to control for pop cases in which there may be insertions/deletions
## or whether this should be controlled by the user
countSites_cherryMethDiff <- function(cherries, data) {
  ## returns a data frame with one line for each cherry, containing the branch distance between the
  ## cherries and for each structure two numbers: that of the number of sites that are m in one in
  ## u in the other tip and the number of sites that are p in one and m or u in the other tip.

  str_n <- length(data[[1]]) # set the number of structures
  df <- data.frame(tips=character(), dist=numeric()) # start df
  for(str in 1:str_n) {
    # e.g. df 3 structures: tips dist 1_f 1_h 2_f 2_h 3_f 3_h
    df <- cbind(df, integer(), integer())
    names(df)[c(ncol(df)-1, ncol(df))] <- paste0(str,"_",c("f", "h"))
  }
  for(i in 1:length(cherries)) {
    ch <- cherries[[i]]
    df[i, 1] <- paste0(ch[1],"-",ch[2]) # add tips info
    df[i, 2] <- ch[3] # add dist between tips
    for(str in 1:str_n) {
      # count number of sites with full methylation change (from u to m or m to u)
      df[i, (str+1)*2-1] <- sum(abs(data[[ch[1]]][[str]]-data[[ch[2]]][[str]])==1)
      # count number of sites with half methylation change (p in one tip and m or u in the other)
      df[i, (str+1)*2] <- sum(abs(data[[ch[1]]][[str]]-data[[ch[2]]][[str]])==0.5)
    }
  }
  df
}

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
