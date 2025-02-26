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

#' Validate and Parse a Phylogenetic Tree
#'
#' This function checks whether the input is a valid phylogenetic tree, either as a character string in Newick format
#' or as an object of class \code{phylo} from the \code{ape} package. If the input is a Newick string, it is parsed into 
#' a \code{phylo} object. The function also ensures that the tree contains at least two tips.
#'
#' @param tree A phylogenetic tree in Newick format (as a character string) or an object of class \code{phylo} from the \code{ape} package.
#'   - If the input is a character string, it must follow the Newick format (e.g., \code{"((tip_1:1,tip_2:1):5,tip_3:6);"}).
#'   - If an object of class \code{phylo} is provided, it should be a valid phylogenetic tree.
#'
#' @return A \code{phylo} object representing the validated and parsed tree.
#'
#' @details 
#' - The function first verifies that the input is either a valid \code{phylo} object or a character string.
#' - If the input is a Newick string, it attempts to parse it into a \code{phylo} object using \code{ape::read.tree()}.
#' - If parsing fails, an informative error message is returned.
#' - The function also checks that the tree contains at least two tips, as a valid phylogenetic tree should have at least one split.
#'
#' @importFrom ape read.tree
validate_tree <- function(tree){
  
  # Check that the input is a character string or phylo object
  if (!is.character(tree) && !inherits(tree, "phylo")) {
    stop(paste("Argument 'tree' must be a character string in newick format (e.g. '((tip_1:1,tip_2:1):5,tip_3:6);').\n",
               "or an object of class 'phylo' from the 'ape' package. For more details, see ?phylo or ?ape."))
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
      stop(paste("Error in ape::read.tree(tree): Invalid 'tree' format: ", conditionMessage(w), "\n",
                 "Argument 'tree' must be a character string in newick format (e.g. '((tip_1:1,tip_2:1):5,tip_3:6);').\n",
                 "or an object of class 'phylo' from the 'ape' package. For more details, see ?phylo or ?ape."))
    }, error = function(e) {
      stop(paste("Error in ape::read.tree(tree): Invalid 'tree' format: ", conditionMessage(e), "\n",
                 "Argument 'tree' must be a character string in newick format (e.g. '((tip_1:1,tip_2:1):5,tip_3:6);').\n",
                 "or an object of class 'phylo' from the 'ape' package. For more details, see ?phylo or ?ape."))
    })
  }
  
  # Check given tree has minium two tips
  if (length(tree$tip.label)<2) stop("The input 'tree' must have a minimum of 2 tips.")
  
  # Return tree as phylo object
  tree
}

#' Get Cherry Pair Distances from a Phylogenetic Tree
#'
#' This function computes the pairwise distances between the tips of a phylogenetic tree 
#' that are part of cherries. A cherry is a pair of leaf nodes (also called tips or terminal nodes) 
#' in a phylogenetic tree that share a direct common ancestor. 
#' In other words, if two leaves are connected to the same internal node and no other leaves 
#' are connected to that internal node, they form a cherry.
#' The distance is calculated as the sum of the branch lengths between the two cherry tips.
#'
#' @param tree A tree in Newick format (as a character string) or an object of class \code{phylo} from the \code{ape} package.
#' If the input is a character string, it must follow the Newick or New Hampshire format (e.g. \code{"((tip_1:1,tip_2:1):5,tip_3:6);"}).
#' If an object of class \code{phylo} is provided, it should represent a valid phylogenetic tree.
#'
#' @param input_control A logical value indicating whether to validate the input tree. 
#' If \code{TRUE} (default), the function checks that the tree is in a valid format and has at least two tips.
#' If \code{FALSE}, the function assumes the tree is already valid and skips the validation step.
#'
#' @return A data frame with five columns:
#'   \item{first_tip_name}{A character string representing the name of the first tip in the cherry.}
#'   \item{second_tip_name}{A character string representing the name of the second tip in the cherry.}
#'   \item{first_tip_index}{An integer representing the index of the first tip in the cherry.}
#'   \item{second_tip_index}{An integer representing the index of the second tip in the cherry.}
#'   \item{dist}{A numeric value representing the sum of the branch lengths between the two tips (i.e., the distance between the cherries).}
#'   
#' @details The function first checks if the input is either a character string in the Newick format or an object of class \code{phylo}, 
#' unless \code{input_control} is set to \code{FALSE}. It then computes the pairwise distances between the tips in the tree and 
#' identifies the sister pairs (cherries). The distance between each cherry is the sum of the branch lengths leading to the sister tips.
#'   
#'   The tips of each cherry are identified by their names and indices. 
#'   The tip indices correspond to (a) the index from left to right on the Newick string, 
#'   (b) the order of the tip label in the \code{phylo_object$tip.label}, and 
#'   (c) the index in the methylation data list (\code{data[[tip]][[structure]]}) as obtained with the function \code{simulate_evolData()} when the given tree has several tips.
#'   
#'   If the tree is provided in Newick format, it will be parsed using the \code{ape::read.tree} function.
#'   
#' @importFrom ape read.tree cophenetic.phylo
#'   
#' @examples
#' # Example of a tree in Newick format
#' 
#' newick_tree <- "((a:1,b:2):5,c:6);"
#' 
#' get_cherryDist(newick_tree)
#' 
#' # Example of using a phylo object from ape
#' 
#' library(ape)
#' tree_phylo <- read.tree(text = "((a:1,b:1):5,c:6);")
#' 
#' get_cherryDist(tree_phylo)
#' 
#' @export
get_cherryDist <- function(tree, input_control = TRUE){
  
  # Check input tree format and minium two tips, get tree in phylo format (ape package)
  if (input_control) tree <- validate_tree(tree)
  
  # Compute the pairwise distances between the tips from a phylogenetic tree
  dist <- ape::cophenetic.phylo(tree)
  # Set a vector to save the cherry tips for which the distance has already been extracted (because dist is symmetric)
  cherry_tips <- c()
  # set df to store the cherry names, tip indices and distances
  cherry_dist <- data.frame(first_tip_name=character(), 
                            second_tip_name=character(),
                            first_tip_index=integer(),
                            second_tip_index=integer(),
                            dist=numeric()) 
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
            # save the tip labels, the distance and the tip indices
            cherry_dist[nrow(cherry_dist)+1,] <- c(tip, tip_a, which(rownames(dist)==tip), which(rownames(dist)==tip_a), dist[tip,tip_a])
            cherry_tips <- c(cherry_tips, tip, tip_a)
          }
        }
      }
    }
  }
  cherry_dist$dist <- as.numeric(cherry_dist$dist)
  cherry_dist$first_tip_index <- as.integer(cherry_dist$first_tip_index)
  cherry_dist$second_tip_index <- as.integer(cherry_dist$second_tip_index)
  cherry_dist
}


#' Validate Structure of Input Data for Cherry Distance Computation
#'
#' This function checks whether the provided input data has the required structure.
#' It ensures that the number of tips is sufficient and that the data structure is consistent across tips and structures.
#'
#' @param cherryDist A data frame containing cherry pair distances, including tip indices (output from \code{get_cherryDist})
#' @param data A nested list representing structured data for each tip, following the format \code{data[[tip]][[structure]]}.
#'
#' @details The function performs several validation steps:
#' - Ensures that the number of tips in \code{data} is at least as large as the highest tip index in \code{cherryDist}.
#' - Checks that all tips contain at least one structure and that the number of structures is consistent across tips.
#' - Verifies that within each structure, all tips have the same number of sites and no zero-length structures.
#'
#' If any of these conditions fail, the function throws an error with a descriptive message.
validate_data <- function(cherryDist, data){
  
  # Get the number of tips in data
  num_tips<- length(data)
  
  # Check the number of tips is equal or greater to the maximum tip index in cherryDist argument
  if(!num_tips>=max(c(cherryDist$first_tip_index, cherryDist$second_tip_index))){
    stop("Argument 'data' with required structure data[[tip]][[structure]] does not have enough tips")
  }
  
  # Ensure the number of structures is > 0 and equal across tips
  lengths_per_tip <- sapply(seq_len(num_tips), function(tip) length(data[[tip]]))
  if (!(all(lengths_per_tip > 0) && length(unique(lengths_per_tip)) == 1)) {
    stop("Input 'data' with required structure data[[tip]][[structure]] has some tips with zero or differing number of structures. Check given 'data' structure.")
  }
  
  # Set the number of structures after checking consistency across tips
  str_n <- length(data[[1]]) 
  
  # Iterate over each structure and check length (number of sites) consistency across tips
  for (structure in seq_len(str_n)) {
    lengths_per_tip <- sapply(seq_len(num_tips), function(tip) length(data[[tip]][[structure]]))
    
    # Ensure all lengths are > 0 and equal across tips
    if (!(all(lengths_per_tip > 0) && length(unique(lengths_per_tip)) == 1)) {
      stop(paste("Error: Input 'data' with required structure data[[tip]][[structure]] has inconsistent lengths across tips or zero length at some tips at structure", structure, "."))
    }
  }
}

#' Count Methylation Differences Between Cherry Pairs
#'
#' This function calculates the number of methylation differences between pairs of cherry tips in a phylogenetic tree.
#' A cherry is a pair of leaf nodes that share a direct common ancestor. The function quantifies full and half methylation
#' differences for each genomic structure (e.g., island/non-island) across all sites.
#'
#' @param cherryDist A data frame containing pairwise distances between the tips of a phylogenetic tree that form cherries.
#'   This should be as the output of \code{get_cherryDist}, and must include the following columns:
#'   \describe{
#'     \item{first_tip_name}{A character string representing the name of the first tip in the cherry.}
#'     \item{second_tip_name}{A character string representing the name of the second tip in the cherry.}
#'     \item{first_tip_index}{An integer representing the index of the first tip in the cherry.}
#'     \item{second_tip_index}{An integer representing the index of the second tip in the cherry.}
#'     \item{dist}{A numeric value representing the sum of the branch lengths between the two tips (i.e., the distance between the cherries).}
#'   }
#'
#' @param data A list containing methylation states at tree tips for each genomic structure (e.g., island/non-island).
#'   The data should be structured as \code{data[[tip]][[structure]]}, where each structure has the same number of sites across tips.
#'   The input data must be prefiltered to ensure CpG sites are represented consistently across different tips.
#'   
#' @param input_control A logical value indicating whether to validate the input data. 
#' If \code{TRUE} (default), the function checks that the data has the required structure.
#' It ensures that the number of tips is sufficient and that the data structure is consistent across tips and structures.
#' If \code{FALSE}, the function assumes the tree is already valid and skips the validation step.
#'
#' @return A data frame with one row per cherry, containing the following columns:
#'   \describe{
#'     \item{tip_names}{A character string representing the names of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{tip_indices}{A character string representing the indices of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{dist}{A numeric value representing the sum of the branch distances between the cherry tips.}
#'     \item{One column for each structure named with the structure number followed by _f}{An integer count of the sites with a full methylation difference (where one tip is methylated and the other is unmethylated) for the given structure.}
#'     \item{One column for each structure named with the structure number followed by _h}{An integer count of the sites with a half methylation difference (where one tip is partially methylated and the other is either fully methylated or unmethylated) for the given structure.}
#'   }
#'
#' @details The function first verifies that \code{cherryDist} contains the required columns and has at least one row.
#'   It also ensures that \code{data} contains a sufficient number of tips and that all structures have the same number of sites.
#'   The function then iterates over each cherry and genomic structure to compute the number of full and half methylation differences
#'   between the two tips of each cherry.
#'
#' @examples
#' # Example data setup
#' 
#' data <- list(
#'   list(c(0, 1, 2, 0), c(1, 1, 0, 2)),
#'   list(c(1, 0, 2, 1), c(0, 1, 2, 2))
#' )
#' 
#' tree <- "(tip1:0.25, tip2:0.25);"
#' 
#' cherryDist <- get_cherryDist(tree)
#' 
#' countSites_cherryMethDiff(cherryDist, data)
#'
#' @export
countSites_cherryMethDiff <- function(cherryDist, data, input_control = TRUE) {
  
  # Check argument cherryDist has all columns
  if(!all(c("first_tip_name", "second_tip_name", "first_tip_index", "second_tip_index", "dist") %in% names(cherryDist))){
    stop("Argument 'cherryDist' misses a required column")
  }
  
  # Check argument cherryDist has at least one row
  if(!nrow(cherryDist)>0) stop("Argument cherryDist has 0 rows.")
  
  # Check input data format and minium number of tips
  if (input_control) validate_data(cherryDist, data)
  
  # Set the number of structures after checking consistency across tips
  str_n <- length(data[[1]]) 
  
  # Initialize a data frame to store the count of methylation differences at each cherry
  df <- data.frame(tip_names=character(), tip_indices=character(), dist=numeric()) 
  
  # For each structure in the data add two integer columns named as
  # [structure_index]_f for full methylation changes and [structure_index]_h for half methylation changes
  for(str in 1:str_n) {
    df <- cbind(df, integer(), integer())
    names(df)[c(ncol(df)-1, ncol(df))] <- paste0(str,"_",c("f", "h"))
  }
  
  # Loop through the cherries
  for(i in 1:nrow(cherryDist)) {
    
    # Extract current cherry info
    ch <- cherryDist[i,]
    
    # Save current cherry tip names, tip indices and distance
    df[i, "tip_names"] <- paste0(ch$first_tip_name,"-",ch$second_tip_name) 
    df[i, "tip_indices"] <- paste0(ch$first_tip_index,"-",ch$second_tip_index) 
    df[i, "dist"] <- ch$dist 
    
    # Loop through the structures
    for(str in 1:str_n) {
      # Get the sequence of methylation states for the given structure at the cherry tips
      data_first_tip <- data[[ch$first_tip_index]][[str]]
      data_second_tip <- data[[ch$second_tip_index]][[str]]
      
      # Count number of sites with full methylation change (from u to m or m to u)
      df[i, paste0(str, "_f")] <- sum(abs(data_first_tip-data_second_tip)==1)
      
      # Count number of sites with half methylation change (p in one tip and m or u in the other)
      df[i, paste0(str, "_h")] <- sum(abs(data_first_tip-data_second_tip)==0.5)
    }
  }
  # Return the data frame
  df
}

#' Compute Methylation Frequency Differences Between Cherry Pairs
#'
#' This function calculates the frequency of methylation differences between pairs of cherry tips in a phylogenetic tree.
#' A cherry is a pair of leaf nodes that share a direct common ancestor. The function quantifies full and half methylation
#' differences for each genomic structure (e.g., island/non-island) across all sites and normalizes these counts by the number
#' of sites per structure to obtain frequencies.
#'
#' @param tree A phylogenetic tree object. The function assumes it follows an appropriate format for downstream processing.
#'
#' @param data A list containing methylation states at tree tips for each genomic structure (e.g., island/non-island).
#'   The data should be structured as `data[[tip]][[structure]]`, where each structure has the same number of sites across tips.
#'   The input data must be prefiltered to ensure CpG sites are represented consistently across different tips.
#'
#' @param input_control A logical value indicating whether to validate the input data.
#'   If `TRUE` (default), the function checks that the data has the required structure.
#'   It ensures that the number of tips is sufficient and that the data structure is consistent across tips and structures.
#'   If `FALSE`, the function assumes the tree is already valid and skips the validation step.
#'
#' @return A data frame with one row per cherry, containing the following columns:
#'   \describe{
#'     \item{tip_names}{A character string representing the names of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{tip_indices}{A character string representing the indices of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{dist}{A numeric value representing the sum of the branch distances between the cherry tips.}
#'     \item{One column for each structure named with the structure number followed by _f}{A numeric value representing the frequency of sites with a full methylation difference (where one tip is methylated and the other is unmethylated) for the given structure.}
#'     \item{One column for each structure named with the structure number followed by _h}{A numeric value representing the frequency of sites with a half methylation difference (where one tip is partially methylated and the other is either fully methylated or unmethylated) for the given structure.}
#'   }
#'
#' @details The function first validates the tree structure and extracts pairwise distances between cherry tips.
#'   It then counts methylation differences using `countSites_cherryMethDiff` and normalizes these counts by the number
#'   of sites per structure to compute frequencies. The resulting data frame provides a per-cherry frequency
#'   of methylation differences (half or full difference) across different genomic structures.
#'
#' @examples
#' # Example data setup
#' 
#' data <- list(
#' list(rep(1,10), rep(0,5), rep(1,8)),
#' list(rep(1,10), rep(0.5,5), rep(0,8)),
#' list(rep(1,10), rep(0.5,5), rep(0,8)),
#' list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
#' 
#' tree <- "((a:1.5,b:1.5):2,(c:2,d:2):1.5);"
#' 
#' freqSites_cherryMethDiff(tree, data)
#'
#' @export
freqSites_cherryMethDiff <- function(tree, data, input_control = TRUE){
  
  # Check input tree format and minium two tips, get tree in phylo format (ape package)
  if (input_control) tree <- validate_tree(tree)
  
  # Get cherry distances aoviding duplicate input control
  cherryDist <- get_cherryDist(tree, input_control = FALSE)
  
  # Check input tree data format and minium number of tips
  if (input_control) validate_data(cherryDist, data)
    
  # Get the count numbers per type (full or half) of methylation change per cherry
  # avoiding duplicate input control
  df <- countSites_cherryMethDiff(cherryDist, data, input_control = FALSE)
  
  # Get the number of structures
  str_n <- length(data[[1]])
  
  # Get the number of sites per structure
  sites_n <- numeric(length = str_n)
  for(str in 1:str_n) sites_n[str] <- length(data[[1]][[str]]) 
  
  # Duplicate them so that there is one-to-one correspondence with the counts in df
  dup_str_n <- rep(sites_n, each = 2)
  
  # Filter the counts (3 first columns correspond to tip_names, tip_indices and dist)
  # to be transformed to frequencies per structure and methylation type
  freqs <- df[,4:ncol(df)]
  
  # Divide each value by the number of sites in the corresponding strucute
  for(cherry in 1:nrow(freqs)) freqs[cherry,] <- freqs[cherry,]/dup_str_n
  
  # Update the frequencies in the dataframe containing tip_names, tip_indices and dist
  df[,4:ncol(df)] <- freqs
  
  df
}


#' Compute Site Frequency of Methylation Changes per Cherry
#'
#' This function calculates the total frequency of methylation differences (both full and half changes)
#' for each genomic structure for each cherry in a phylogenetic tree. 
#' A cherry is a pair of leaf nodes (also called tips or terminal nodes) 
#' in a phylogenetic tree that share a direct common ancestor. 
#' In other words, if two leaves are connected to the same internal node and no other leaves 
#' are connected to that internal node, they form a cherry.
#'
#' @param tree A phylogenetic tree in Newick format or a phylo object from the ape package. The function ensures
#'   the tree has a valid structure and at least two tips.
#'
#' @param data A list containing methylation states at tree tips for each genomic structure
#'   (e.g., island/non-island). The data should be structured as \code{data[[tip]][[structure]]}, where
#'   each structure has the same number of sites across tips.
#'
#' @return A data frame with one row per cherry, containing the following columns:
#'   \describe{
#'     \item{tip_names}{A character string representing the names of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{tip_indices}{A character string representing the indices of the two tips in the cherry, concatenated with a hyphen.}
#'     \item{dist}{A numeric value representing the sum of the branch distances between the cherry tips.}
#'     \item{One column for each structure named with the structure number}{A numeric value representing the total frequency of methylation changes (both full and half) for the given structure.}
#'   }
#'
#' @details The function first verifies that \code{tree} and \code{data} have valid structures and the minimum number of tips.
#'   It then extracts per-cherry methylation differences using \code{freqSites_cherryMethDiff}, handling potential errors.
#'   Finally, it aggregates the full and half methylation differences for each genomic structure at each cherry.
#'
#' @examples
#' # Example data setup
#' 
#' data <- list(
#' list(rep(1,10), rep(0,5), rep(1,8)),
#' list(rep(1,10), rep(0.5,5), rep(0,8)),
#' list(rep(1,10), rep(0.5,5), rep(0,8)),
#' list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
#' 
#' tree <- "((a:1.5,b:1.5):2,(c:2,d:2):1.5);"
#' 
#' get_siteFChange_cherry(tree, data)
#'
#' @export
get_siteFChange_cherry <- function(tree, data){
  
  # Check input tree format and minium two tips, get tree in phylo format (ape package),
  # check input tree data format and minium number of tips, and 
  # get per-cherry frequency of (half and full) methylation differences across different genomic structures
  tryCatch({
    freqSites_perMethDiffType <- freqSites_cherryMethDiff(tree, data)
  }, warning = function(w) {
    stop(conditionMessage(w))
  }, error = function(e) {
    stop(conditionMessage(e))
  })

  # Identify the unique structures in the column names (3 first columns correspond to tip_names, tip_indices and dist)
  structures <- unique(sub("_[fh]$", "", names(freqSites_perMethDiffType)[-c(1, 2, 3)]))
  
  # Create the new dataframe
  siteFChange_cherry <- data.frame(
    tip_names = freqSites_perMethDiffType$tip_names,
    tip_indices = freqSites_perMethDiffType$tip_indices,
    dist = freqSites_perMethDiffType$dist
  )
  
  # Sum for each structure the total frequency of methylation change (both half and full)
  for (str in structures){
    siteFChange_cherry[[str]] <- rowSums(freqSites_perMethDiffType[, grep(paste0("^", str, "_"), names(freqSites_perMethDiffType))])
  }
  
  siteFChange_cherry
}








