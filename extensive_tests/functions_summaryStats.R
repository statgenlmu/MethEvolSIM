



#### #### #### Fraction of methylated sites over methylated and unmethylated #### #### ####

##### Mean fraction of methylated sites over methylated and unmethylated sites in islands
## index_islands: vector with structural indeces for islands
## data: list with methylation states at tree tips for each structure 
## data[[tip]][[structure]] when the number of tips is >1, 
## or data[[structure]] when there is only one tip. Methylation states are vectors
## sample_n: number of given tips/samples
get_islandMeanFracMoverMU <- function(index_islands, data, sample_n){
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
      mean_tip[s] <- mean(data[[s]][[i]]==1)/(1-mean(data[[s]][[i]]==0.5))
    }
    mean_island[island_counter] <- mean(mean_tip, na.rm = TRUE)
    island_counter <- island_counter + 1
  }
  return(mean(mean_island, na.rm = TRUE))
}

##### Mean fraction of methylated sites over methylated and unmethylated sites in non-islands
## index_nonislands: vector with structural indeces for non-islands
## data: list with methylation states at tree tips for each structure 
## data[[tip]][[structure]] when the number of tips is >1, 
## or data[[structure]] when there is only one tip. Methylation states are vectors
## sample_n: number of given tips/samples
get_nonislandMeanFracMoverMU <- function(index_nonislands, data, sample_n){
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
      mean_tip[s] <- mean(data[[s]][[i]]==1)/(1-mean(data[[s]][[i]]==0.5))
    }
    mean_nonisland[nonisland_counter] <- mean(mean_tip, na.rm = TRUE)
    nonisland_counter <- nonisland_counter + 1
  }
  return(mean(mean_nonisland, na.rm = TRUE))
}

##### Standard deviation fraction of methylated sites over methylated and unmethylated sites in islands
## index_islands: vector with structural indeces for islands
## data: list with methylation states at tree tips for each structure 
## data[[tip]][[structure]] when the number of tips is >1, 
## or data[[structure]] when there is only one tip. Methylation states are vectors
## sample_n: number of given tips/samples
get_islandSDFracMoverMU <- function(index_islands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    frac_island <- c()
    island_counter <- 1
    for (i in index_islands){
      # Compute proportion of partially methylated sites at each island
      frac_island[island_counter] <- mean(data[[s]][[i]] == 1)/(1-mean(data[[s]][[i]]==0.5))
      island_counter <- island_counter + 1
    }
    sd_tip[s] <- sd(frac_island, na.rm = TRUE)
  }
  return(mean(sd_tip))
}


##### Standard deviation fraction of methylated sites over methylated and unmethylated sites in non-islands
## index_nonislands: vector with structural indeces for non-islands
## data: list with methylation states at tree tips for each structure 
## data[[tip]][[structure]] when the number of tips is >1, 
## or data[[structure]] when there is only one tip. Methylation states are vectors
## sample_n: number of given tips/samples
get_nonislandSDFracMoverMU <- function(index_nonislands, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  sd_tip <- c() 
  for (s in 1:sample_n){
    frac_nonisland <- c()
    nonisland_counter <- 1
    for (i in index_nonislands){
      # Compute proportion of partially methylated sites at each island
      frac_nonisland[nonisland_counter] <- mean(data[[s]][[i]] == 1)/(1-mean(data[[s]][[i]]==0.5))
      nonisland_counter <- nonisland_counter + 1
    }
    sd_tip[s] <- sd(frac_nonisland, na.rm = TRUE)
  }
  return(mean(sd_tip))
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

cherry_diffs <- function(cherries, data) {
  ## returns a data frame with one line for each cherry, containing the branch distance between the
  ## cherries and for each structure two numbers: that of the number of sites that are m in one in
  ## u in the other tip and the number of sites that are p in one and m or u in the other tip.
  n <- length(data)[[1]] # set the number of structures
  df <- data.frame(tips=character(), dist=numeric())
  for(j in 1:n) {
    df <- cbind(df, integer(), integer())
    names(df)[c(ncol(df)-1, ncol(df))] <- paste0(j,"_",c("f", "h"))
  }
  for(i in 1:length(cherries)) {
    ch <- cherries[[i]]
    df[i, 1] <- paste0(ch[1],"-",ch[2])
    df[i, 2] <- ch[3]
    for(j in 1:n) {
      df[i, (j+1)*2-1] <- sum(abs(data[[ch[1]]][[j]]-data[[ch[2]]][[j]])==1)
      df[i, (j+1)*2] <- sum(abs(data[[ch[1]]][[j]]-data[[ch[2]]][[j]])==0.5)
    }
  }
  df
}

get_FChange_cherryData <- function(tree, data, n_CpG){
  cherries <- find_cherries(tree)
  FChange_cherryData <- cherry_diffs(cherries = cherries, data = data)
  # Each value is a count number of full methylation changes _f or changes in one strand _h.
  # Each structure has a number of CpGs. Then, if we wanna count the proportion
  # of sites, we divide each value by the number of CpGs. 
  # Then, as each _f case means a change in the 2 strands and each _h case in one,
  # we multiply each _f value by 2.
  FChange_cherryData[,c(3:ncol(FChange_cherryData))] <- FChange_cherryData[,c(3:ncol(FChange_cherryData))]/n_CpG
  FChange_cherryData[,grep("f", colnames(FChange_cherryData))] <- FChange_cherryData[,grep("f", colnames(FChange_cherryData))]*2
  return(FChange_cherryData)
}

get_meanFChange_cherry_i <- function(FChange_cherryData, index_islands){
  FChange_cherryData$meanFChange_cherry_i <- rowMeans(FChange_cherryData[,c(paste0(index_islands, "_f"), paste0(index_islands, "_h"))])
  return(FChange_cherryData[,c("dist", "meanFChange_cherry_i")])
}

get_meanFChange_cherry_ni <- function(FChange_cherryData, index_nonislands){
  FChange_cherryData$meanFChange_cherry_ni <- rowMeans(FChange_cherryData[,c(paste0(index_nonislands, "_f"), paste0(index_nonislands, "_h"))])
  return(FChange_cherryData[,c("dist", "meanFChange_cherry_ni")])
}

fit_meanFChange_cherry_i <- function(FChange_cherryData, index_islands){
  FChange_cherryData_iMean <- get_meanFChange_cherry_i(FChange_cherryData, index_islands)
  start_a <- max(FChange_cherryData_iMean$meanFChange_cherry_i)
  mod <- nls(meanFChange_cherry_i ~ a * (1 - exp(-r *dist)), 
             data = FChange_cherryData_iMean,
             start = list(a = start_a, r = 0.1),
             control = nls.control(maxiter = 1000))
  mod_s <- summary(mod)
  return(list(rate = mod_s$coefficients["r",1],
              scale = mod_s$coefficients["a",1]))
}

fit_meanFChange_cherry_ni <- function(FChange_cherryData, index_nonislands){
  FChange_cherryData_niMean <- get_meanFChange_cherry_ni(FChange_cherryData, index_nonislands)
  start_a <- max(FChange_cherryData_niMean$meanFChange_cherry_ni)
  mod <- nls(meanFChange_cherry_ni ~ a * (1 - exp(-r *dist)), 
             data = FChange_cherryData_niMean,
             start = list(a = start_a, r = 0.1),
             control = nls.control(maxiter = 1000))
  mod_s <- summary(mod)
  return(list(rate = mod_s$coefficients["r",1],
              scale = mod_s$coefficients["a",1]))
}


## TODO: Add argument to script -- development
development <- TRUE
if(development){
  library(devtools)
  load_all()
} else {
  library(MethEvolSIM) ## to steal function MethEvolSIM:::split_newick
}


## tree: tree in newick format
## meth: array methylation data at tips of the tree, vectors of u, p, m, rows indexed with tip names as in tree
## returns list of [[1]] a list of sets of optimal states for the root of tree for each island (or position) and
##    [[2]] a vector of minimum number of changes needed for each island (or position)
compute_fitch <- function(tree, meth) {
  if(substr(tree,1,1)!="(") {
    ## seems to be a tip
    return(list(meth[as.integer(tree), ], rep(0, dim(meth)[2])))
  }
  ## if still here, then real tree
  subtree <- MethEvolSIM:::split_newick(tree)$unit
  if(length(subtree)!=2) stop("function compute_fitch works only for rooted binary trees")
  L <- compute_fitch(subtree[1], meth)
  R <- compute_fitch(subtree[2], meth)
  minchan <- L[[2]] + R[[2]]
  optset <- list()
  for(i in 1:length(minchan)) {
    optset[[i]] <- intersect(L[[1]][[i]], R[[1]][[i]])
    if(length(optset[[i]])==0) {
      optset[[i]] <- union(L[[1]][[i]], R[[1]][[i]])
      minchan[i] <- minchan[i] + 1
    } 
  }
  return(list(optset, minchan))
}

categorize <- function(x, uthr, mthr) {
  ## convert methylation frequencies of uthr or lower as "u", of "mthr" or higher as "m" and in between as "p"
  s <- character()
  for(i in 1:length(x)) {
    if(x[[i]] <= uthr) {
      s[i] <- "u"
    } else {
      if(x[[i]] >= mthr) {
        s[i] <- "m"
      } else {
        s[i] <- "p"
      }
    }
  }
  return(s)
}

computeFitch_RegionGlbSt <- function(index_islands, data, tree, u_threshold, m_threshold){
  ## average methylation of CpG islands:
  cpgavg <- lapply(data, function(x) sapply(x[index_islands], mean))
  ## convert averaged methylation to categories "u", "p", and "m"
  upmdata_pre <- lapply(cpgavg, categorize, u_threshold, m_threshold)
  ## Reorder tips according to the tip labels in the tree;
  upmdata <- upmdata_pre[order(as.integer(tree$tip.label))]
  upmdata  <- matrix(unlist(upmdata), nrow=100, byrow=TRUE)
  result <- compute_fitch(ape::write.tree(tree), upmdata)[[2]]
}

#### #### #### Transitions between DMR #### #### ####


##### Select the borders between regions with a difference in mean methylation levels over a threshold
select_minDiffMethTrans <- function(threshold, data, sample_n) {
  
  # Restructure data as nested list when data corresponds to a single replicate/tip/sample
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }

  minDiffMethTrans <- data.frame(border = integer(),
                                 tip = integer(),
                                 type = factor(levels = c(0, 1), labels = c("i", "d")))
  # Set the number of borders
  border_n <- length(data[[1]]) - 1
  #print(paste("border_n", border_n))
  for (border in 1:border_n) {
    diff_mean <- numeric(length(data)) # vector to store differences in mean structural methylation
    
    for (tip in 1:length(data)) {
      left_mean <- mean(data[[tip]][[border]])
      right_mean <- mean(data[[tip]][[border + 1]])
      diff_mean[tip] <- left_mean - right_mean
    }
    
    # Find the tips where the absolute difference is greater than the threshold
    index_minDiff <- which(abs(diff_mean) > threshold)
    
    for (selected_tip in index_minDiff) {
      type <- ifelse(diff_mean[selected_tip] > 0, "d", "i")
      
      # Append to data frame
      minDiffMethTrans <- rbind(minDiffMethTrans, data.frame(border = border,
                                                             tip = selected_tip,
                                                             type = factor(type, levels = c("i", "d"))))
    }
  }
  
  return(minDiffMethTrans)
}

select_minRepresentation <- function(minDiffMethTrans, minRepresentation){
  selected_data <- list(border = c(),
                        tips = list())
  # Create a table of counts of border type per border
  type_counts <- table(minDiffMethTrans$border, minDiffMethTrans$type)
  # Create a list to store the average methylation at the borders ## TODO:
  counter <- 1
  for(border in unique(minDiffMethTrans$border)){
    # Determine the most common type for the current border and tips with that transition type
    border_counts <- type_counts[as.character(border), ]
    most_common_type <- names(which.max(border_counts))
    selected_tips <- minDiffMethTrans$tip[minDiffMethTrans$border == border & minDiffMethTrans$type== most_common_type]
    # Average the transition methylations for the given number of CpGs at each side of the transiton
    if (length(selected_tips) >= minRepresentation){
      selected_data$border[counter] <- border
      selected_data$tips[[counter]]<- selected_tips
      counter <- counter + 1
    }
  }
  return(selected_data)
} 


## TODO: delete this one
compute_kTrans <- function(border, tips, subset_CpG_n){
  # Set vectors sum the methylation states at each side of the border for the selected tips
  sum_leftStr <- rep(0, subset_CpG_n)
  sum_rightStr <- rep(0, subset_CpG_n)
  for(tip in tips){
    sum_leftStr <- sum_leftStr + data[[tip]][[border]][(length(data[[tip]][[border]])-(subset_CpG_n-1)):length(data[[tip]][[border]])]
    sum_rightStr <- sum_rightStr + data[[tip]][[border+1]][1:subset_CpG_n]
  }
  # Transform methylation values (0,0.5 or 1) to count values (double strand: 0, 1 or 2)
  k <- 2*c(sum_leftStr, sum_rightStr)
  return(k)
}

## TODO: document this one
compute_kTrans_acrossSamples <- function(data, border, tips, subset_CpG_n){
  
  # Set vectors sum the methylation states at each side of the border for the selected tips
  sum_leftStr <- rep(0, subset_CpG_n)
  sum_rightStr <- rep(0, subset_CpG_n)
  for(tip in tips){
    minCpG_check <- check_minCpGnumber(data, border, subset_CpG_n, sample_n = length(data))
    if(minCpG_check){
      for(tip in tips){
        sum_leftStr <- sum_leftStr + data[[tip]][[border]][(length(data[[tip]][[border]])-(subset_CpG_n-1)):length(data[[tip]][[border]])]
        sum_rightStr <- sum_rightStr + data[[tip]][[border+1]][1:subset_CpG_n]
      }
    }
  }
  # Transform methylation values (0,0.5 or 1) to count values (double strand: 0, 1 or 2)
  k <- 2*c(sum_leftStr, sum_rightStr)
  return(k)
}

## TODO: document this one
compute_kTrans_withinSample <- function(data, minDiffMethTrans, subset_CpG_n){
  # Set vectors sum the methylation states at each side of the border for the selected tips
  sum_leftStr <- rep(0, subset_CpG_n)
  sum_rightStr <- rep(0, subset_CpG_n)
  for (i in 1:nrow(minDiffMethTrans)){
    border <- minDiffMethTrans$border[i]
    minCpG_check <- check_minCpGnumber(data, border, subset_CpG_n, sample_n = 1)
    if(minCpG_check){
      if(minDiffMethTrans$type[i] == "d"){ # Reverse the methylation states at borders with decreasing transitions
        sum_leftStr <- sum_leftStr + rev(data[[border]][(length(data[[border]])-(subset_CpG_n-1)):length(data[[border]])])
        sum_rightStr <- sum_rightStr + rev(data[[border+1]][1:subset_CpG_n])
      } else {
        sum_leftStr <- sum_leftStr + data[[border]][(length(data[[border]])-(subset_CpG_n-1)):length(data[[border]])]
        sum_rightStr <- sum_rightStr + data[[border+1]][1:subset_CpG_n]
      }
    }
  }
  # Transform methylation values (0,0.5 or 1) to count values (double strand: 0, 1 or 2)
  k <- 2*c(sum_leftStr, sum_rightStr)
  return(k)
}

## TODO: document this one
check_minCpGnumber <- function(data, border, subset_CpG_n, sample_n){
  # Restructure data as nested list when data corresponds to a single replicate/tip/samle
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  if(length(data[[1]][[border]])<2*subset_CpG_n || length(data[[1]][[border+1]])<2*subset_CpG_n){
    print(paste("Structures at the side of each border should have at least twice the number of CpGs than given 'subset_CpG_n'. Neglecting border:", border))
    return(FALSE)
  } else {
    TRUE
  }
}


# inverse composite likelihood function to optimize # heuristic approximation
minusCLL <- function(params, x, k, n) {
  left_limit <- params[1]
  right_limit <- params[2]
  expected_midpoint <- params[3]
  steepness <- params[4]
  
  # Calculate the logistic function values
  logistic <- left_limit + (right_limit - left_limit) / (1 + exp(-steepness * (x - expected_midpoint)))
  
  # Check for non-finite logistic values
  if (any(!is.finite(logistic))) {
    print(paste("Non-finite logistic values for params:", paste(params, collapse=", ")))
    return(Inf)
  }
  
  # Calculate the log-likelihood
  # Removed the log(n choose k) because it does not depend on the logistic
  CLL <- sum(k * log(logistic) + (n - k) * log(1 - logistic))
  
  # Return negative log-likelihood
  return(-CLL)
}

fit_logistic <- function(k, n, subset_CpG_n){
  # adjust the initial values for the parameters to be optimized
  left_limit <- mean(k[1:(subset_CpG_n/2)]/n)
  right_limit <- mean(k[(length(k)-(subset_CpG_n/2)-1):length(k)]/n)
  expected_midpoint <- subset_CpG_n + 0.5
  steepness <- 1.5 
  param <- c(left_limit, right_limit, expected_midpoint, steepness)
  # Set the x values to fit the function
  x <- 1:(2*subset_CpG_n)
  fitted_params <- optim(par=param, minusCLL, lower=c(0.01, 0.01, 10, 0.1), upper=c(0.99, 0.99, 50, 10),
                         x=x, k=k, n=n, method="L-BFGS-B")$par
  return(list(left_limit = fitted_params[1],
              right_limit = fitted_params[2],
              midpoint = fitted_params[3],
              steepness = fitted_params[4]))
}

## TODO: Update comment.
##### Get the mean and sd of the steepness of fitted logistic transition between DMR
## data: list with methylation states at tree tips for each structure 
## data[[tip]][[structure]] when the number of tips is >1
## threshold: minimum difference of average methylation frequencies between regions to select the border. Needs to be smaller than the number of tips.
## minRepresentation: minimum number of samples/replicates with a selected border of equal type (decreasing or increasing methylation)
## subset_CpG_n: number of CpG sites at each side of the border to consider.
fit_MethTrans <- function(data, threshold, minRepresentation, subset_CpG_n, sample_n){
  minDiffMethTrans <- select_minDiffMethTrans(threshold, data, sample_n)
  steepness <- c()
  midpoint <- c()
  if (sample_n == 1){
    # Use all the borders of minDiffMethTrans if there is a minimum number of borders (minRepresentation).
    # For that, reverse decreasing transitions (type == "d")
    if(minRepresentation > length(data)-1){
      stop("The minimum number of borders to select transitions 'minRepresentation' is higher than the number of borders in the data")
    }
    if(nrow(minDiffMethTrans) >= minRepresentation){
      # fit
      k <- compute_kTrans_withinSample(data, minDiffMethTrans, subset_CpG_n)
      n <- 2*nrow(minDiffMethTrans) 
      fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
      steepness[1] <- fitted_params$steepness
      midpoint[1] <- fitted_params$midpoint
    } else {
      print("minRepresentation not met")
    }
    
  } else {
    if(!is.list(data[[1]])){
      stop("To fit the methylation transitions data across tips/samples/replicates data needs to be a nested list (e.g: data[[tip]][[structure]])")
    }
    if(minRepresentation > length(data)){
      stop("The minimum number of tips/samples/replicates set to select transitions 'minRepresentation' is higher than the number of tips/samples/replicates in the data")
    }
    # For each border, select the transition type (decreasing or incresing) more represented across samples
    # and if it is represented in a minimum number of replicates (minRepresentation), select it
    selectedTrans <- select_minRepresentation(minDiffMethTrans = minDiffMethTrans, minRepresentation = 10)
    for (i in 1:length(selectedTrans$border)){
      border <- selectedTrans$border[i]
      k <- compute_kTrans_acrossSamples(data, border, tips = selectedTrans$tips[[i]], subset_CpG_n)
      n <- 2*length(selectedTrans$tips[[i]])
      fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
      steepness[i] <- fitted_params$steepness
      midpoint[i] <- fitted_params$midpoint
    }
  }
  ## TODO: also return mean and sd of midpoint
  return(list(meanSteepness = mean(steepness),
              sdSteepness = sd(steepness),
              meanMidPoint = mean(midpoint),
              sdMidPoint = sd(midpoint)))
}

## TODO: After visualizations of summary statistics show new functions work, delete the ones inside if (FALSE)
if (FALSE) {
  ## TODO: Update comment.
  ##### Get the mean and sd of the steepness of fitted logistic transition between DMR
  ## data: list with methylation states at tree tips for each structure 
  ## data[[tip]][[structure]] when the number of tips is >1
  ## threshold: minimum difference of average methylation frequencies between regions to select the border. Needs to be smaller than the number of tips.
  ## minRepresentation: minimum number of samples/replicates with a selected border of equal type (decreasing or increasing methylation)
  ## subset_CpG_n: number of CpG sites at each side of the border to consider.
  fit_MethTrans <- function(data, threshold, minRepresentation, subset_CpG_n, sample_n){
    minDiffMethTrans <- select_minDiffMethTrans(threshold, data, sample_n)
    steepness <- c()
    midpoint <- c()
    if (sample_n == 1){
      # Use all the borders of minDiffMethTrans if there is a minimum number of borders (minRepresentation).
      # For that, reverse decreasing transitions (type == "d")
      if(minRepresentation > length(data)-1){
        stop("The minimum number of borders to select transitions 'minRepresentation' is higher than the number of borders in the data")
      }
      if(nrow(minDiffMethTrans) >= minRepresentation){
        # fit
        k <- compute_kTrans_withinSample(data, minDiffMethTrans, subset_CpG_n)
        n <- 2*nrow(minDiffMethTrans) 
        fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
        steepness[i] <- fitted_params$steepness
        midpoint[i] <- fitted_params$midpoint
      } else {
        print("minRepresentation not met")
      }
      
    } else {
      if(!is.list(data[[1]])){
        stop("To fit the methylation transitions data across tips/samples/replicates data needs to be a nested list (e.g: data[[tip]][[structure]])")
      }
      if(minRepresentation > length(data[[1]])){
        stop("The minimum number of tips/samples/replicates set to select transitions 'minRepresentation' is higher than the number of tips/samples/replicates in the data")
      }
      # For each border, select the transition type (decreasing or incresing) more represented across samples
      # and if it is represented in a minimum number of replicates (minRepresentation), select it
      selectedTrans <- select_minRepresentation(minDiffMethTrans = minDiffMethTrans, minRepresentation = 10)
      for (i in 1:length(selectedTrans$border)){
        k <- compute_kTrans_acrossSamples(data, border, tips = selectedTrans$tips[[i]], subset_CpG_n)
        n <- 2*length(selectedTrans$tips[[i]])
        fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
        steepness[i] <- fitted_params$steepness
        midpoint[i] <- fitted_params$midpoint
      }
    }
    ## TODO: also return mean and sd of midpoint
    return(list(meanSteepness = mean(steepness),
                sdSteepness = sd(steepness)))
  }
  
  ## TODO: This one wil be fit_MethTrans_acrossSamples
  ##### Get the mean and sd of the steepness of fitted logistic transition between DMR
  ## data: list with methylation states at tree tips for each structure 
  ## data[[tip]][[structure]] when the number of tips is >1
  ## threshold: minimum difference of average methylation frequencies between regions to select the border. Needs to be smaller than the number of tips.
  ## minRepresentation: minimum number of samples/replicates with a selected border of equal type (decreasing or increasing methylation)
  ## subset_CpG_n: number of CpG sites at each side of the border to consider.
  fit_MethTrans <- function(data, threshold, minRepresentation, subset_CpG_n){
    if(!is.list(data[[1]])){
      stop("To fit the methylation transitions data at several tips/samples/replicates needs to be provided (e.g: data[[tip]][[structure]])")
    }
    if(minRepresentation > length(data[[1]])){
      stop("The minimum number of tips/samples/replicates set to select transitions 'minRepresentation' is higher than the number of tips/samples/replicates in the data")
    }
    minDiffMethTrans <- select_minDiffMethTrans(threshold, data)
    selectedTrans <- select_minRepresentation(minDiffMethTrans = minDiffMethTrans, minRepresentation = minRepresentation)
    steepness <- c()
    midpoint <- c()
    for (i in 1:length(selectedTrans$border)){
      border <- selectedTrans$border[i]
      if(length(data[[1]][[border]])<2*subset_CpG_n || length(data[[1]][[border+1]])<2*subset_CpG_n){
        print(paste("Structures at the side of each border should have at least twice the number of CpGs than given 'subset_CpG_n'. Neglecting border:", border))
      } else {
        k <- compute_kTrans(border = border, tips = selectedTrans$tips[[i]], subset_CpG_n = subset_CpG_n)
        n <- 2*length(selectedTrans$tips[[i]])
        fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
        steepness[i] <- fitted_params$steepness
        midpoint[i] <- fitted_params$midpoint
      }
    }
    return(list(meanSteepness = mean(steepness),
                sdSteepness = sd(steepness)))
  }
  
  ## TODO: This one wil be fit_MethTrans_withinSample
  ##### Get the mean and sd of the steepness of fitted logistic transition between DMR
  ## data: list with methylation states at tree tips for each structure 
  ## data[[tip]][[structure]] when the number of tips is >1
  ## threshold: minimum difference of average methylation frequencies between regions to select the border. Needs to be smaller than the number of tips.
  ## minRepresentation: minimum number of samples/replicates with a selected border of equal type (decreasing or increasing methylation)
  ## subset_CpG_n: number of CpG sites at each side of the border to consider.
  fit_MethTrans <- function(data, threshold, minRepresentation, subset_CpG_n){
    if(!is.list(data[[1]])){
      stop("To fit the methylation transitions data at several tips/samples/replicates needs to be provided (e.g: data[[tip]][[structure]])")
    }
    if(minRepresentation > length(data[[1]])){
      stop("The minimum number of tips/samples/replicates set to select transitions 'minRepresentation' is higher than the number of tips/samples/replicates in the data")
    }
    minDiffMethTrans <- select_minDiffMethTrans(threshold, data)
    selectedTrans <- select_minRepresentation(minDiffMethTrans = minDiffMethTrans, minRepresentation = minRepresentation)
    steepness <- c()
    midpoint <- c()
    for (i in 1:length(selectedTrans$border)){
      border <- selectedTrans$border[i]
      if(length(data[[1]][[border]])<2*subset_CpG_n || length(data[[1]][[border+1]])<2*subset_CpG_n){
        print(paste("Structures at the side of each border should have at least twice the number of CpGs than given 'subset_CpG_n'. Neglecting border:", border))
      } else {
        k <- compute_kTrans(border = border, tips = selectedTrans$tips[[i]], subset_CpG_n = subset_CpG_n)
        n <- 2*length(selectedTrans$tips[[i]])
        fitted_params <- fit_logistic(k = k, n = n, subset_CpG_n = subset_CpG_n)
        steepness[i] <- fitted_params$steepness
        midpoint[i] <- fitted_params$midpoint
      }
    }
    return(list(meanSteepness = mean(steepness),
                sdSteepness = sd(steepness)))
  }
}







### OLD EXPLORED SUMMARY STATISTICS. FUNCTIONS NOT UPDATED ###


compute_meanCov_i <- function(index_islands, subset_CpG_n, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  str_counter <- 1
  cov_i <- c()
  for (tip in 1:sample_n){
    for(i in index_islands){
      cov_i[str_counter] <- cov(data[[tip]][[i]][(length(data[[tip]][[i]])-subset_CpG_n):(subset_CpG_n-1)], data[[tip]][[i]][(length(data[[tip]][[i]])-subset_CpG_n +1):(subset_CpG_n)])
      str_counter <- str_counter + 1
    }
  }
  return(mean(cov_i, na.rm = TRUE))
}


compute_meanCov_ni <- function(index_nonislands, subset_CpG_n, data, sample_n){
  # Restructure data as nested list
  if (sample_n == 1){
    data_list <- list()
    data_list[[1]] <- data
    data <- data_list
  }
  str_counter <- 1
  cov_ni <- c()
  for (tip in 1:sample_n){
    for(i in index_nonislands){
      cov_ni[str_counter] <- cov(data[[tip]][[i]][(length(data[[tip]][[i]])-subset_CpG_n):(subset_CpG_n-1)], data[[tip]][[i]][(length(data[[tip]][[i]])-subset_CpG_n +1):(subset_CpG_n)])
      str_counter <- str_counter + 1
    }
  }
  return(mean(cov_ni, na.rm = TRUE))
}


#### #### #### Tip pairwise comparisons #### #### ####
library(ape)

##### Get a dataframe with distances between each combination of pairs of tips
## tree: phylogenetic tree object (as generated with package scrm)
get_pairwiseDistance <- function(tree){
  
  # Get the number of tree tips
  n_tips <- length(tree$tip.label)
  
  # Generate the vectors for the tip pairwise comparisons
  tipA <- numeric()
  tipB <- numeric()
  
  for (tA in 1:(n_tips - 1)){
    for (tB in (tA+1):n_tips){
      tipA <- c(tipA, tA)
      tipB <- c(tipB, tB)
    }
  }
  
  # Get the matrix with the distances between the tree tips
  distance <- signif(ape::cophenetic.phylo(tree), 4)
  
  # Set a dataframe for the pairwise distances
  pairwiseD <- data.frame(tipA = tipA,
                          tipB = tipB,
                          distance = rep(NA, length(tipA)))
  
  # Save the distances
  for (tip in 1:length(tipA)){
    pairwiseD$distance[tip] <- distance[tipA[tip],tipB[tip]]
  }
  pairwiseD
}

##### Compute the mean proportion of sites with a different methylation state in islands between two given tree tips
## index_islands: vector with structural indeces for islands
## tipA: numerical index for the first tip
## tipB: numerical index for the second tip
## data: list with methylation states at tree tips for each structure data[[tip]][[structure]]. Methylation states are vectors
get_meanFChangeSt_i <- function(index_islands, tipA, tipB, data){
  FChangeSt_i <- c()
  str_counter <- 1
  for (i in index_islands){
    # Get the pair of methylation states
    str_tipA <- data[[tipA]][[i]]
    str_tipB <- data[[tipB]][[i]]
    # Calculate the proportion of sites with a different methylation state
    FChangeSt_i[str_counter] <- mean(str_tipA != str_tipB)
    str_counter <- str_counter + 1
  }
  return(mean(FChangeSt_i))
}

##### Compute the mean proportion of sites with a different methylation state in non-islands between two given tree tips
## index_nonislands: vector with structural indeces for non-islands
## tipA: numerical index for the first tip
## tipB: numerical index for the second tip
## data: list with methylation states at tree tips for each structure data[[tip]][[structure]]. Methylation states are vectors
get_meanFChangeSt_ni <- function(index_nonislands, tipA, tipB, data){
  FChangeSt_i <- c()
  str_counter <- 1
  for (i in index_nonislands){
    # Get the pair of methylation states
    str_tipA <- data[[tipA]][[i]]
    str_tipB <- data[[tipB]][[i]]
    # Calculate the proportion of sites with a different methylation state
    FChangeSt_i[str_counter] <- mean(str_tipA != str_tipB)
    str_counter <- str_counter + 1
  }
  return(mean(FChangeSt_i))
}

##### Compute the proportion of tip pairwise comparisons between islands with a difference in proportion of any of the methylation states over a threshold
## index_islands: vector with structural indeces for islands
## tipA: numerical index for the first tip
## tipB: numerical index for the second tip
## upper_threshold: Upper quantile (from 0 to 100)
## data: list with methylation states at tree tips for each structure data[[tip]][[structure]]. Methylation states are vectors
get_diffFreqSt_i <- function(index_islands, tipA, tipB, upper_threshold, data){
  diffFStU_i <- c()
  diffFStP_i <- c()
  diffFStM_i <- c()
  str_counter <- 1
  for (i in index_islands){
    # Get the pair of methylation states
    str_tipA <- data[[tipA]][[i]]
    str_tipB <- data[[tipB]][[i]]
    # Calculate the frequency of the each methylation state for the island at each tip
    freqStU_tipA <- sum(str_tipA == 0)/length(str_tipA)
    freqStU_tipB <- sum(str_tipB == 0)/length(str_tipB)
    freqStP_tipA <- sum(str_tipA == 0.5)/length(str_tipA)
    freqStP_tipB <- sum(str_tipB == 0.5)/length(str_tipB)
    freqStM_tipA <- sum(str_tipA == 1)/length(str_tipA)
    freqStM_tipB <- sum(str_tipB == 1)/length(str_tipB)
    # Calculate the absolute difference in frequency between the tips
    diffFStU_i[str_counter] <- abs((freqStU_tipA - freqStU_tipB)*100)
    diffFStP_i[str_counter] <- abs((freqStP_tipA - freqStP_tipB)*100)
    diffFStM_i[str_counter] <- abs((freqStM_tipA - freqStM_tipB)*100)
    str_counter <- str_counter + 1
  }
  
  # Calculate the proportion of times the absolute difference in frequency is equal or higher than the threshold
  prop_overT <- (mean(diffFStU_i >= upper_threshold) + mean(diffFStP_i >= upper_threshold) + mean(diffFStM_i >= upper_threshold))/3
  return(prop_overT)
}

##### Compute the proportion of tip pairwise comparisons between non-islands with a difference in proportion of any of the methylation states over a threshold
## index_nonislands: vector with structural indeces for non-islands
## tipA: numerical index for the first tip
## tipB: numerical index for the second tip
## upper_threshold: Upper quantile (from 0 to 100)
## data: list with methylation states at tree tips for each structure data[[tip]][[structure]]. Methylation states are vectors
get_diffFreqSt_ni <- function(index_nonislands, tipA, tipB, upper_threshold, data){
  diffFStU_ni <- c()
  diffFStP_ni <- c()
  diffFStM_ni <- c()
  str_counter <- 1
  for (i in index_nonislands){
    # Get the pair of methylation states
    str_tipA <- data[[tipA]][[i]]
    str_tipB <- data[[tipB]][[i]]
    # Calculate the frequency of the each methylation state for the island at each tip
    freqStU_tipA <- sum(str_tipA == 0)/length(str_tipA)
    freqStU_tipB <- sum(str_tipB == 0)/length(str_tipB)
    freqStP_tipA <- sum(str_tipA == 0.5)/length(str_tipA)
    freqStP_tipB <- sum(str_tipB == 0.5)/length(str_tipB)
    freqStM_tipA <- sum(str_tipA == 1)/length(str_tipA)
    freqStM_tipB <- sum(str_tipB == 1)/length(str_tipB)
    # Calculate the absolute difference in frequency between the tips
    diffFStU_ni[str_counter] <- abs((freqStU_tipA - freqStU_tipB)*100)
    diffFStP_ni[str_counter] <- abs((freqStP_tipA - freqStP_tipB)*100)
    diffFStM_ni[str_counter] <- abs((freqStM_tipA - freqStM_tipB)*100)
    str_counter <- str_counter + 1
  }
  
  # Calculate the proportion of times the absolute difference in frequency is equal or higher than the threshold
  prop_overT <- (mean(diffFStU_ni >= upper_threshold) + mean(diffFStP_ni >= upper_threshold) + mean(diffFStM_ni >= upper_threshold))/3
  return(prop_overT)
}



compare_tips <- function(tree, index_islands, index_nonislands, data){
  pairwiseComp <- get_pairwiseDistance(tree)
  pairwiseComp$meanFChangeSt_i <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$meanFChangeSt_ni <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_1_i <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_1_ni <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_5_i <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_5_ni <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_50_i <- rep(NA, nrow(pairwiseComp))
  pairwiseComp$diffFSt_50_ni <- rep(NA, nrow(pairwiseComp))
  
  for(comp in 1:nrow(pairwiseComp)){
    pairwiseComp$meanFChangeSt_i[comp] <- get_meanFChangeSt_i(index_islands = index_islands, 
                                                              tipA = pairwiseComp$tipA[comp], 
                                                              tipB = pairwiseComp$tipB[comp], 
                                                              data = data)
    pairwiseComp$meanFChangeSt_ni[comp] <- get_meanFChangeSt_ni(index_nonislands = index_nonislands, 
                                                                tipA = pairwiseComp$tipA[comp], 
                                                                tipB = pairwiseComp$tipB[comp], 
                                                                data = data)
    pairwiseComp$diffFSt_1_i[comp] <- get_diffFreqSt_i(index_islands = index_islands, 
                                                       tipA = pairwiseComp$tipA[comp], 
                                                       tipB = pairwiseComp$tipB[comp], 
                                                       upper_threshold = 1,
                                                       data = data)
    pairwiseComp$diffFSt_1_ni[comp] <- get_diffFreqSt_ni(index_nonislands = index_nonislands, 
                                                         tipA = pairwiseComp$tipA[comp], 
                                                         tipB = pairwiseComp$tipB[comp], 
                                                         upper_threshold = 1,
                                                         data = data)
    pairwiseComp$diffFSt_5_i[comp] <- get_diffFreqSt_i(index_islands = index_islands, 
                                                       tipA = pairwiseComp$tipA[comp], 
                                                       tipB = pairwiseComp$tipB[comp], 
                                                       upper_threshold = 5,
                                                       data = data)
    pairwiseComp$diffFSt_5_ni[comp] <- get_diffFreqSt_ni(index_nonislands = index_nonislands, 
                                                         tipA = pairwiseComp$tipA[comp], 
                                                         tipB = pairwiseComp$tipB[comp], 
                                                         upper_threshold = 5,
                                                         data = data)
    pairwiseComp$diffFSt_50_i[comp] <- get_diffFreqSt_i(index_islands = index_islands, 
                                                        tipA = pairwiseComp$tipA[comp], 
                                                        tipB = pairwiseComp$tipB[comp], 
                                                        upper_threshold = 50,
                                                        data = data)
    pairwiseComp$diffFSt_50_ni[comp] <- get_diffFreqSt_ni(index_nonislands = index_nonislands, 
                                                          tipA = pairwiseComp$tipA[comp], 
                                                          tipB = pairwiseComp$tipB[comp], 
                                                          upper_threshold = 50,
                                                          data = data)
  }
  return(pairwiseComp)
}




