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


