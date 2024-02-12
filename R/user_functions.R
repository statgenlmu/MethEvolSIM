# Island equilibrium frequencies sampling strategia
singleStructureGenerator$set("public", "get_alpha_pI", function() private$alpha_pI)
singleStructureGenerator$set("public", "get_beta_pI", function() private$beta_pI)
singleStructureGenerator$set("public", "get_alpha_mI", function() private$alpha_mI)
singleStructureGenerator$set("public", "get_beta_mI", function() private$beta_mI)
# Non-island equilibrium frequencies sampling strategia
singleStructureGenerator$set("public", "get_alpha_pNI", function() private$alpha_pNI)
singleStructureGenerator$set("public", "get_beta_pNI", function() private$beta_pNI)
singleStructureGenerator$set("public", "get_alpha_mNI", function() private$alpha_mNI)
singleStructureGenerator$set("public", "get_beta_mNI", function() private$beta_mNI)
# DNA methylation state evolution
## IWE
combiStructureGenerator$set("public", "get_mu", function() private$mu)
## SSE
singleStructureGenerator$set("public", "get_alpha_Ri", function() private$alpha_Ri)
singleStructureGenerator$set("public", "get_iota", function() private$iota)


#' Get Default Parameter Values
#'
#' This function retrieves parameter values for the DNA methylation simulation.
#'
#' @param rootData NULL to return default parameter values. For data parameter values, provide rootData$data as the output of simulate_initialData().
#'
#' @return A data frame containing default parameter values.
#'
#' @details The function called without arguments returns default parameter values.
#' When rootData (as $data output of simulate_initialData()) is given, it returns data parameter values.
#'
#' @export
get_parameterValues <- function(rootData = NULL){
  if(is.null(rootData)){
    obj <- singleStructureGenerator$new("U", 10)
    infoStr <- data.frame(start = c(1, 14, 15),
                          end = c(13, 14, 19),
                          globalState = c("M", "M", "M"))
    combi_obj <- combiStructureGenerator$new(infoStr)
  } else {
    if(class(rootData)[1] != "combiStructureGenerator"){
      stop("'rootData' needs to be the $data output of simulate$initialData")
    } else {
      obj <- rootData$get_singleStr(1)
      combi_obj <- rootData
    }
  }
  data.frame(alpha_pI = obj$get_alpha_pI(),
             beta_pI = obj$get_beta_pI(),
             alpha_mI = obj$get_alpha_mI(),
             beta_mI = obj$get_beta_mI(),
             alpha_pNI = obj$get_alpha_pNI(),
             beta_pNI = obj$get_beta_pNI(),
             alpha_mNI = obj$get_alpha_mNI(),
             beta_mNI = obj$get_beta_mNI(),
             mu = combi_obj$get_mu(),
             alpha_Ri = obj$get_alpha_Ri(),
             iota = obj$get_iota())
}

#' Simulate Initial Data
#'
#' This function simulates initial data based on the provided information and parameters.
#'
#' @param infoStr A data frame containing columns 'start', 'end', and 'globalState'.
#'  If customized equilibrium frequencies are given, it also contains columns 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq'
#'  with the equilibrium frequency values for unmethylated, partially methylated and methylated.
#' @param params Optional data frame with specific parameter values.
#' Structure as in get_parameterValues() output. If not provided, default values will be used.
#'
#' @return A list containing the simulated data ($data) and parameters ($params).
#'
#' @details The function performs several checks on the input data and parameters
#'  to ensure they meet the required criteria and simulates DNA methylation data.
#'
#' @export
simulate_initialData <- function(infoStr, params = NULL){
  if (!is.data.frame(infoStr) ||
      !all(c("start", "end", "globalState") %in% colnames(infoStr))) {
    stop("infoStr should be a dataframe with columns: 'start', 'end', 'globalState'")
  }
  if(all(c("u_eqFreq", "p_eqFreq", "m_eqFreq") %in% colnames(infoStr))){
    for (i in 1:nrow(infoStr)){
      eqFreqs <- c(infoStr$u_eqFreq[i], infoStr$p_eqFreq[i], infoStr$m_eqFreq[i])
      if(any(is.na(eqFreqs))){
        stop(paste("if 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq' are given, they need to be frequency values. Missing values in row ", i))
      } else if(!is.numeric(eqFreqs) || !length(eqFreqs) == 3 || !sum(eqFreqs)==1){
        stop(paste("if 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq' are given, they need to be frequency values. Incorrect values in row ", i))
      }
    }
  }
  if(!is.null(params)){
    if(!is.data.frame(params) || !all(c("alpha_pI", "beta_pI", "alpha_mI", "beta_mI", "alpha_pNI", "beta_pNI", "alpha_mNI", "beta_mNI", "mu", "alpha_Ri", "iota") %in% colnames(params))){
      stop("if 'params' is given, it needs to be a dataframe with column names as in get_parameterValues() output")
    }
  }
  data <- combiStructureGenerator$new(infoStr = infoStr, params = params)
  if(is.null(params)){
    params <- get_parameterValues()
  }
  return(list(data = data,
              params = params))
}

#' Simulate Data Evolution along a Tree
#'
#' This function simulates DNA methylation data either by simulating data at the root of the provided evolutionary tree
#' (if infoStr is given) or by using pre-existing data at the root (if rootData is given) and letting it evolve along the tree.
#'
#' @param infoStr A data frame containing columns 'start', 'end', and 'globalState'.
#'  If customized initial equilibrium frequencies are given, it also contains columns 'u_eqFreq', 'p_eqFreq', and 'm_eqFreq'
#'  with the equilibrium frequency values for unmethylated, partially methylated, and methylated.
#' @param rootData The output of the simulate_initialData() function. It represents the initial data at the root of the evolutionary tree.
#' @param tree A string in Newick format representing the evolutionary tree.
#' @param params Optional data frame with specific parameter values.
#' Structure as in get_parameterValues() output. If not provided, default values will be used.
#' @param dt Length of time step for the simulation (default is 0.01).
#'
#' @return A list containing the simulated data ($data), parameters ($params), the evolutionary tree ($tree), and the time step for the simulation ($dt).
#'
#' @details The function performs several checks on the input data and parameters
#'  to ensure they meet the required criteria and simulates DNA methylation data along the provided evolutionary tree.
#'
#' @export
simulate_evolData <- function(infoStr = NULL, rootData = NULL, tree, params = NULL, dt = 0.01){
  if (is.null(infoStr) && is.null(rootData)) {
    stop("At least one of infoStr or rootData must be provided.")
  } else if (!is.null(infoStr) && !is.null(rootData)) {
    stop("Only one of infoStr or rootData should be provided.")
  }
  if(!is.null(rootData)){
    if(!is.null(params)){
      stop("When rootData is given, rootData parameter values are used. Argument 'params' needs to be null. \n To get rootData parameter values use get_parameterValues(rootData). \n To customize rootData parameter values use simulate_initialData()")
    }
    print("Parameter values set as in given rootData")
    if(class(rootData)[1] != "combiStructureGenerator"){
      stop("rootData should be the output of simulate_initialData()")
    }
  }
  if(!is.null(infoStr)){
    if (!is.data.frame(infoStr) ||
        !all(c("start", "end", "globalState") %in% colnames(infoStr))) {
      stop("infoStr should be a dataframe with columns: 'start', 'end', 'globalState'")
    }
  }
  if(!is.character(tree)){
    stop("tree needs to be given as string in newick format")
  }
  if(all(c("u_eqFreq", "p_eqFreq", "m_eqFreq") %in% colnames(infoStr))){
    for (i in 1:nrow(infoStr)){
      eqFreqs <- c(infoStr$u_eqFreq[i], infoStr$p_eqFreq[i], infoStr$m_eqFreq[i])
      if(any(is.na(eqFreqs))){
        stop(paste("if 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq' are given, they need to be frequency values. Missing values in row ", i))
      } else if(!is.numeric(eqFreqs) || !length(eqFreqs) == 3 || !sum(eqFreqs)==1){
        stop(paste("if 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq' are given, they need to be frequency values. Incorrect values in row ", i))
      }
    }
  }
  if(!is.null(params)){
    if(!is.data.frame(params) || !all(c("alpha_pI", "beta_pI", "alpha_mI", "beta_mI", "alpha_pNI", "beta_pNI", "alpha_mNI", "beta_mNI", "mu", "alpha_Ri", "iota") %in% colnames(params))){
      stop("if 'params' is given, it needs to be a dataframe with column names as in get_parameterValues() output")
    }
  }
  data <- treeMultiRegionSimulator$new(infoStr = infoStr, rootData = rootData, tree = tree, params = params)
  if(is.null(params)){
    params <- get_parameterValues(rootData)
  }
  return(list(data = data,
              params = params,
              tree = tree,
              dt = dt))
}


