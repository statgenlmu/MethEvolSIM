#### #### #### Frequency of methylation states per structure type #### #### ####

#' Calculate the Mean Frequency of Partially Methylated Sites in Islands
#'
#' This function computes the mean frequency of partially methylated sites (with methylation state 0.5) 
#' for the set of genomic structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each genomic structure (island / non-island) 
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean frequency of partially methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
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
#' for a set of genomic structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each genomic structure (island / non-island).
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
#' index_islands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
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
#' (with methylation state 0.5) for a set of genomic structures identified as islands.
#'
#' @param index_islands A vector containing the structural indices for islands.
#' @param data A list containing methylation states at tree tips for each genomic structure (island / non-island).
#' For a single tip: \code{data[[structure]]}. For multiple tips: \code{data[[tip]][[structure]]}.
#' @param sample_n The number of samples (tips) to process.
#'
#' @return A numeric value representing the mean standard deviation of partially methylated sites in the islands.
#' @examples
#' # Example usage:
#' index_islands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
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
#' (with methylation state 0.5) for a set of genomic structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each genomic structure (island / non-island).
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
#' (with methylation state 1) for a set of genomic structures identified as islands.
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
#' index_islands <- c(1, 3)
#' data <- list(
#'   list(c(0.5, 1, 0.5), c(0, 0.5, 1), c(1, 0, 0.5)), # tip 1
#'   list(c(0.5, 0.5, 0), c(1, 0.5, 0.5), c(0.5, 0.5, 1)) # tip 2
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
#' (with methylation state 1) for a set of genomic structures identified as non-islands.
#'
#' @param index_nonislands A vector containing the structural indices for non-islands.
#' @param data A list containing methylation states at tree tips for each genomic structure (island / non-island).
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


#### #### #### Mean neighbor correlations per structure type #### #### ####



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

#' Get Cherry Pair Distances from a Phylogenetic Tree
#'
#' This function computes the pairwise distances between the tips of a phylogenetic tree 
#' that are part of cherries. A cherry is a pair of leaves (tips) that are sisters in the tree.
#' The distance is calculated as the sum of the branch lengths between the two sister tips.
#'
#' @param tree A tree in Newick format (as a character string) or an object of class \code{phylo} from the \code{ape} package.
#' If the input is a character string, it must follow the Newick or New Hampshire format (e.g. \code{"((tip_1:1,tip_2:1):5,tip_3:6);"}).
#' If an object of class \code{phylo} is provided, it should represent a valid phylogenetic tree.
#'
#' @return A data frame with three columns:
#'   \item{first_tip}{A character string representing the first tip in the cherry (sister pair).}
#'   \item{second_tip}{A character string representing the second tip in the cherry.}
#'   \item{dist}{A numeric value representing the sum of the branch lengths between the two tips (i.e., the distance between the cherries).}
#'   
#' @details The function first checks if the input is either a character string in the Newick format or an object of class \code{phylo}.
#'   It then computes the pairwise distances between the tips in the tree and identifies the sister pairs (cherries).
#'   The distance between each cherry is the sum of the branch lengths leading to the sister tips.
#'   
#'   If the tree is provided in Newick format, it will be parsed using the \code{ape::read.tree} function.
#'   
#' @importFrom ape read.tree cophenetic.phylo
#'   
#' @examples
#' # Example of a tree in Newick format
#' newick_tree <- "((a:1,b:2):5,c:6);"
#' get_cherryDist(newick_tree)
#' 
#' # Example of using a phylo object from ape
#' library(ape)
#' tree_phylo <- read.tree(text = "((a:1,b:1):5,c:6);")
#' get_cherryDist(tree_phylo)
#' 
#' @export
get_cherryDist <- function(tree){
  
  # Check that the input is a character string or phylo object
  if (!is.character(tree) && !inherits(tree, "phylo")) {
    stop(paste("Argument 'tree' must be a character string in newick format (e.g. '((tip_1:1,tip_2:1):5,tip_3:6);').\n",
               "or an object of class 'phylo' from the 'ape'. For more details, see ?phylo or ?ape."))
  }
  
  # Transform the newick tree in a phylo object managing issues with given tree format
  if(!inherits(tree, "phylo")) {
    tryCatch({
      # Attempt to read the tree
      tree <- ape::read.tree(text = tree)
      # Check if the result is NULL, which indicates an issue with the tree format
      if (is.null(tree)) {
        stop("Error in ape::read.tree(tree): Tree could not be parsed. Invalid format.")
      }
    }, warning = function(w) {
      stop("Error in ape::read.tree(tree): Invalid 'tree' format: ", conditionMessage(w))
    }, error = function(e) {
      stop("Error in ape::read.tree(tree): Invalid 'tree' format: ", conditionMessage(e))
    })
  }
  
  # Check given tree has minium two tips
  if (length(tree$tip.label)<2) stop("The input 'tree' must have a minimum of 2 tips.")
  
  # Compute the pairwise distances between the tips from a phylogenetic tree
  dist <- ape::cophenetic.phylo(tree)
  # Set a vector to save the cherry tips for which the distance has already been extracted (because dist is symmetric)
  cherry_tips <- c()
  # set list to store the cherry distances
  cherry_dist <- data.frame(first_tip=character(), second_tip=character(), dist=numeric()) # start df
  # Get the tip names
  tips <- rownames(dist)

  for(tip in tips) {
    # If tip is not already in cherry_dist
    if (!tip %in% cherry_tips) {
      # find the closest tip
      tip_a <- names(which(dist[tip,] == min(dist[tip, colnames(dist) != tip])))
      # if there is one single closest tip
      if(length(tip_a)==1){
        # find the closes tip 
        tip_b <- names(which(dist[tip_a,] == min(dist[tip_a, colnames(dist) != tip_a])))
        # if there is one single closest tip
        if (length(tip_b)==1){
          # if both are mutually closest tips 
          if(tip_b == tip){
            # save the distance and the tip labels 
            cherry_dist[nrow(cherry_dist)+1,] <- c(tip, tip_a, dist[tip,tip_a])
            cherry_tips <- c(cherry_tips, tip, tip_a)
          }
        }
      }
    }
  }
  cherry_dist$dist <- as.numeric(cherry_dist$dist)
  cherry_dist
}