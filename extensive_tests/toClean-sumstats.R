#### #### #### Frequency of methylation states per structure type #### #### ####

#' Calculate the Mean Frequency of Partially Methylated Sites in Islands
#'
#' This function computes the mean frequency of partially methylated sites (with methylation state 0.5) 
#' for a set of structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island) 
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean frequency of partially methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 2)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5)) # tip 2
#' )
#' sample_n <- 2
#' get_islandMeanFreqP(index_islands, data, sample_n)
#' 
#' @export
get_islandMeanFreqP <- function(index_islands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  mean_island <- c()
  island_counter <- 1
  for (i in index_islands){
    mean_tip <- c()
    for(s in 1:sample_n){
      mean_tip[s] <- mean(data[[s]][[i]]==0.5)
    }
    mean_island[island_counter] <- mean(mean_tip)
    island_counter <- island_counter + 1
  }
  return(mean(mean_island))
}

#' Calculate the Mean Frequency of Partially Methylated Sites in Non-Islands
#'
#' This function computes the mean frequency of partially methylated sites (with methylation state 0.5) 
#' for a set of structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean frequency of partially methylated sites in the non-islands.
#' @examples
#' # Example usage:
#' index_nonislands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
#' )
#' sample_n <- 2
#' get_nonislandMeanFreqP(index_nonislands, data, sample_n)
#' 
#' @export
get_nonislandMeanFreqP <- function(index_nonislands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  mean_nonisland <- c()
  nonisland_counter <- 1
  for (i in index_nonislands){
    mean_tip <- c()
    for(s in 1:sample_n){
      mean_tip[s] <- mean(data[[s]][[i]]==0.5)
    }
    mean_nonisland[nonisland_counter] <- mean(mean_tip)
    nonisland_counter <- nonisland_counter + 1
  }
  return(mean(mean_nonisland))
}

#' Calculate the Mean Frequency of Methylated Sites in Islands
#'
#' This function computes the mean frequency of methylated sites (with methylation state 1) 
#' for a set of structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island). 
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean frequency of methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 2)
#' data <- list(
#'   list(c(1, 0, 1), c(1, 1, 0.5)), # tip 1
#'   list(c(1, 0.5, 1), c(0, 1, 1)) # tip 2
#' )
#' sample_n <- 2
#' get_islandMeanFreqM(index_islands, data, sample_n)
#' 
#' @export
get_islandMeanFreqM <- function(index_islands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  mean_island <- c()
  island_counter <- 1
  for (i in index_islands){
    mean_tip <- c()
    for(s in 1:sample_n){
      mean_tip[s] <- mean(data[[s]][[i]]==1)
    }
    mean_island[island_counter] <- mean(mean_tip)
    island_counter <- island_counter + 1
  }
  return(mean(mean_island))
}

#' Calculate the Mean Frequency of Methylated Sites in Non-Islands
#'
#' This function computes the mean frequency of methylated sites (with methylation state 1) 
#' for a set of structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island). 
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean frequency of methylated sites in the non-islands.
#' @examples
#' # Example usage:
#' index_nonislands <- c(1, 3)
#' data <- list(
#'   list(c(1, 0, 1), c(0.5, 1, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(1, 0.5, 1), c(0.5, 1, 1), c(1, 0.5, 0.5)) # tip 2
#' )
#' sample_n <- 2
#' get_nonislandMeanFreqM(index_nonislands, data, sample_n)
#' 
#' @export
get_nonislandMeanFreqM <- function(index_nonislands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  mean_nonisland <- c()
  nonisland_counter <- 1
  for (i in index_nonislands){
    mean_tip <- c()
    for(s in 1:sample_n){
      mean_tip[s] <- mean(data[[s]][[i]]==1)
    }
    mean_nonisland[nonisland_counter] <- mean(mean_tip)
    nonisland_counter <- nonisland_counter + 1
  }
  return(mean(mean_nonisland))
}

#' Calculate the Mean Standard Deviation of Partially Methylated Sites in Islands
#'
#' This function computes the mean standard deviation of partially methylated sites 
#' (with methylation state 0.5) for a set of structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean standard deviation of partially methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 2)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5)) # tip 2
#' )
#' sample_n <- 2
#' get_islandSDFreqP(index_islands, data, sample_n)
#' 
#' @export
get_islandSDFreqP <- function(index_islands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    freq_island <- c()
    island_counter <- 1
    for (i in index_islands){
      # Compute proportion of partially methylated sites at each island
      freq_island[island_counter] <- mean(data[[s]][[i]] == 0.5)
      island_counter <- island_counter + 1
    }
    sd_tip[s] <- stats::sd(freq_island)
  }
  return(mean(sd_tip))
}

#' Calculate the Mean Standard Deviation of Partially Methylated Sites in Non-Islands
#'
#' This function computes the mean standard deviation of partially methylated sites 
#' (with methylation state 0.5) for a set of structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean standard deviation of partially methylated sites in the non-islands.
#' @examples
#' # Example usage:
#' index_nonislands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
#' )
#' sample_n <- 2
#' get_nonislandSDFreqP(index_nonislands, data, sample_n)
#' 
#' @export
get_nonislandSDFreqP <- function(index_nonislands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    freq_nonisland <- c()
    nonisland_counter <- 1
    for (i in index_nonislands){
      # Compute proportion of partially methylated sites at each island
      freq_nonisland[nonisland_counter] <- mean(data[[s]][[i]] == 0.5)
      nonisland_counter <- nonisland_counter + 1
    }
    sd_tip[s] <- stats::sd(freq_nonisland)
  }
  return(mean(sd_tip))
}

#' Calculate the Mean Standard Deviation of Methylated Sites in Islands
#'
#' This function computes the mean standard deviation of methylated sites 
#' (with methylation state 1) for a set of structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' Methylation states are represented as vectors.
#' @param sample_n The number of tips (samples) to process.
#'
#' @return A numeric value representing the mean standard deviation of methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 2)
#' data <- list(
#'   list(c(1, 0.5, 1), c(0, 1, 1)), # tip 1
#'   list(c(1, 1, 0), c(1, 1, 0.5)) # tip 2
#' )
#' sample_n <- 2
#' get_islandSDFreqM(index_islands, data, sample_n)
#' 
#' @export
get_islandSDFreqM <- function(index_islands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    freq_island <- c()
    island_counter <- 1
    for (i in index_islands){
      # Compute proportion of methylated sites at each island
      freq_island[island_counter] <- mean(data[[s]][[i]] == 1)
      island_counter <- island_counter + 1
    }
    sd_tip[s] <- stats::sd(freq_island)
  }
  return(mean(sd_tip))
}


#' Calculate the Mean Standard Deviation of Methylated Sites in Non-Islands
#'
#' This function computes the mean standard deviation of methylated sites 
#' (with methylation state 1) for a set of structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' Methylation states are represented as vectors.
#' @param sample_n The number of tips (samples) to process.
#'
#' @return A numeric value representing the mean standard deviation of methylated sites in the non-islands.
#' @examples
#' # Example usage:
#' index_nonislands <- c(1, 3)
#' data <- list(
#'   list(c(1, 1, 1), c(0, 1, 0.5), c(1, 0, 1)), # tip 1
#'   list(c(1, 0.5, 0), c(1, 1, 0.5), c(1, 1, 1)) # tip 2
#' )
#' sample_n <- 2
#' get_nonislandSDFreqM(index_nonislands, data, sample_n)
#' 
#' @export
get_nonislandSDFreqM <- function(index_nonislands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    freq_nonisland <- c()
    nonisland_counter <- 1
    for (i in index_nonislands){
      # Compute proportion of methylated sites at each non-island
      freq_nonisland[nonisland_counter] <- mean(data[[s]][[i]] == 1)
      nonisland_counter <- nonisland_counter + 1
    }
    sd_tip[s] <- stats::sd(freq_nonisland)
  }
  return(mean(sd_tip))
}



#' Compute the Mean Correlation of Methylation State in Islands
#'
#' This function calculates the mean correlation of methylation states within 
#' island structures, allowing to exclude the shores.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param minN_CpG The minimum number of central CpGs required for computation.
#' @param shore_length The number of CpGs at each side of an island to exclude (shores).
#' @param data A list containing methylation states at tree tips for each structure.
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' Methylation states are represented as vectors.
#' @param sample_n The number of tips (samples) to process.
#'
#' @return A numeric value representing the mean correlation of methylation states in the central CpGs of islands.
#' @details The function processes only islands with a minimum length equal to \code{2 * shore_length + minN_CpG}. 
#' If none has minimum length, returns NA
#' @examples
#' # Example usage:
#' index_islands <- c(1, 2)
#' data <- list(
#'   list(c(0, 1, 0.5, 1, 0.5, 0), c(0.5, 0.5, 1, 1, 0, 0)), # tip 1
#'   list(c(1, 0, 1, 1, 0.5, 0), c(1, 1, 0.5, 0.5, 0, 1))   # tip 2
#' )
#' minN_CpG <- 2
#' shore_length <- 1
#' sample_n <- 2
#' compute_meanCor_i(index_islands, minN_CpG, shore_length, data, sample_n)
#' 
#' @export
compute_meanCor_i <- function(index_islands, minN_CpG, shore_length, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  str_counter <- 1
  cor <- c(NA) 
  for (tip in 1:sample_n){
    for (i in index_islands){
      if(length(data[[tip]][[i]]) >= 2*shore_length + minN_CpG){
        # Define start and end indices to extract the middle segment
        start1 <- shore_length + 1
        start2 <- start1 + 1
        end2 <- length(data[[tip]][[i]]) - shore_length
        end1 <- end2 - 1
        # Extract the sequence info
        segment1 <- data[[tip]][[i]][start1:end1]
        segment2 <- data[[tip]][[i]][start2:end2]
        # Compute the correlation of those segments with methylation state variation
        if (stats::sd(segment1) != 0 && stats::sd(segment2) != 0) {
          cor[str_counter] <- cor(segment1, segment2)
          str_counter <- str_counter + 1
        }
      }
      
    }
  }
  return(mean(cor))
}

#### #### #### Mean correlation per structure type #### #### ####

#' Compute the Mean Correlation of Methylation State in Non-islands
#'
#' This function calculates the mean correlation of methylation states within 
#' non-island structures, allowing to exclude the shores.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param minN_CpG The minimum number of central CpGs required for computation.
#' @param shore_length The number of CpGs at each side of an non-island to exclude (shores).
#' @param data A list containing methylation states at tree tips for each structure.
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' Methylation states are represented as vectors.
#' @param sample_n The number of tips (samples) to process.
#'
#' @return A numeric value representing the mean correlation of methylation states in the central CpGs of non-islands.
#' @details The function processes only non-islands with a minimum length equal to \code{2 * shore_length + minN_CpG}. 
#' If none has minimum length, returns NA
#' @examples
#' # Example usage:
#' index_nonislands <- c(1, 2)
#' data <- list(
#'   list(c(0, 1, 0.5, 1, 0.5, 0), c(0.5, 0.5, 1, 1, 0, 0)), # tip 1
#'   list(c(1, 0, 1, 1, 0.5, 0), c(1, 1, 0.5, 0.5, 0, 1))   # tip 2
#' )
#' minN_CpG <- 2
#' shore_length <- 1
#' sample_n <- 2
#' compute_meanCor_i(index_nonislands, minN_CpG, shore_length, data, sample_n)
#' 
#' @export
compute_meanCor_ni <- function(index_nonislands, minN_CpG, shore_length, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  str_counter <- 1
  cor <- c(NA) 
  for (tip in 1:sample_n){
    for (i in index_nonislands){
      if(length(data[[tip]][[i]]) >= 2*shore_length + minN_CpG){
        # Define start and end indices to extract the middle segment
        start1 <- shore_length + 1
        start2 <- start1 + 1
        end2 <- length(data[[tip]][[i]]) - shore_length
        end1 <- end2 - 1
        # Extract the sequence info
        segment1 <- data[[tip]][[i]][start1:end1]
        segment2 <- data[[tip]][[i]][start2:end2]
        # Compute the correlation of those segments with methylation state variation
        if (stats::sd(segment1) != 0 && stats::sd(segment2) != 0) {
          cor[str_counter] <- cor(segment1, segment2)
          str_counter <- str_counter + 1
        }
      }
      
    }
  }
  return(mean(cor))
}

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
