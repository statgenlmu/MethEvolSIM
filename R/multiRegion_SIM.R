#' @title singleStructureGenerator
#' @importFrom R6 R6Class
#'
#' @description
#' an R6 class representing a single genomic structure
singleStructureGenerator <-
  R6::R6Class("singleStructureGenerator",
              portable = FALSE,
              private = list(
                ## @field my_combiStructure Private attribute: Object of class combiStructureGenerator containing it
                my_combiStructure = NULL,
                ## @field combiStructure_index Private attribute: Index in Object of class combiStructureGenerator
                combiStructure_index = NULL,
                ## @field alpha_pI Private attribute: Model parameter for sampling island equilibrium frequencies
                alpha_pI = 0.1,
                ## @field beta_pI Private attribute: Model parameter for sampling island equilibrium frequencies
                beta_pI = 1,
                ## @field alpha_mI Private attribute: Model parameter for sampling island equilibrium frequencies
                alpha_mI = 0.1,
                ## @field beta_mI Private attribute: Model parameter for sampling island equilibrium frequencies
                beta_mI = 0.5,
                ## @field alpha_pNI Private attribute: Model parameter for sampling non-island equilibrium frequencies
                alpha_pNI = 0.1,
                ## @field beta_pNI Private attribute: Model parameter for sampling non-island equilibrium frequencies
                beta_pNI = 1,
                ## @field alpha_mNI Private attribute: Model parameter for sampling non-island equilibrium frequencies
                alpha_mNI = 0.5,
                ## @field beta_mNI Private attribute: Model parameter for sampling non-island equilibrium frequencies
                beta_mNI = 0.1,
                ## @field globalState Private attribute: Structure's favored global state: "M" for methylated (island structures) / "U" for unmethylated (non-island structures)
                globalState = NULL,
                ## @field eqFreqs Private attribute: Structure's methylation state equilibrium frequencies (for unmethylated, partially methylated and methylated)
                eqFreqs = NULL,

                ## @description
                ## Private method: Sample equilibrium frequencies
                ##
                ## This function samples the object's eqFreqs based on the object's globalState.
                ## Cases when the sample_eqFreqs method is called:
                ##
                ## 1. **During Initialization:** If equilibrium frequencies (eqFreqs) are not provided during the initialization of the object, this method is automatically called to set eqFreqs based on the object's globalState.
                ##
                ## 2. **In IWE_evol():** The method is called within the IWE_evol() function to sample new equilibrium frequencies for the structure.
                ##
                ## @return A numeric vector representing equilibrium frequencies for unmethylated, partially methylated, and methylated states.
                ##
                sample_eqFreqs = function() {
                  if (private$globalState == "M") {
                    p <- stats::rbeta(1, private$alpha_pNI, private$beta_pNI)
                    m <- stats::rbeta(1, private$alpha_mNI, private$beta_mNI) * (1 - p)
                    u <- 1 - p - m; if (u == 1) p = m = 0
                  } else if (private$globalState == "U") {
                    p <- stats::rbeta(1, private$alpha_pI, private$beta_pI)
                    m <- stats::rbeta(1, private$alpha_mI, private$beta_mI) * (1 - p)
                    u <- 1 - p - m; if (u == 1) p = m = 0
                  } else {
                    stop("Invalid globalState")
                  }
                  return(c(u, p, m))
                },
                ## @field seq Private attribute: Encoded sequence of CpGs methylation state: 1 for unmethylated, 2 for partially-methylated, 3 for methylated
                seq = NULL,
                ## @field siteR Private attribute: Encoded sequence of CpGs site SSEi rate: 1 for Ri_values[1], 2 for Ri_values[2] and 3 for Ri_values[3]
                siteR = NULL,
                ## @field neighbSt Private attribute: Encoded sequence of CpGs site neighbor state: as in mapNeighbSt_matrix
                neighbSt = NULL,
                ## @field mapNeighbSt_matrix Private attribute: Matrix encoding neighbor state
                ## Map from neighbor state to numerical coding #########################
                #### #### Rows represent left neighbor state (u,p,m)
                #### #### columns represent right neighbor state (u,p,m)
                #### ####     [,1] [,2] [,3]
                #### ####[1,]    1    2    3
                #### ####[2,]    4    5    6
                #### ####[3,]    7    8    9
                mapNeighbSt_matrix = matrix(c(1L:9L), byrow = TRUE, nrow = 3),
                ## @description
                ## Private method: Get left singleStructureGenerator object in my_combiStructure object last seq state
                ## @return last seq value
                get_leftStr_neighbSt = function(){
                  if (private$combiStructure_index == 1) {
                    stop("'combiStructure_index' is 1, no leftStr available")
                  }
                  if (private$combiStructure_index > 1){
                    private$my_combiStructure$get_singleStr(private$combiStructure_index - 1)$get_seqLastPos()
                  }
                },
                ## @description
                ## Private method: Get right singleStructureGenerator object in my_combiStructure object first seq state
                ## @return first seq value
                get_rightStr_neighbSt = function(){
                  if (private$combiStructure_index == private$my_combiStructure$get_singleStr_number()) {
                    stop("'combiStructure_index' is last, no rightStr available")
                  }
                  if (private$combiStructure_index < private$my_combiStructure$get_singleStr_number()){
                    private$my_combiStructure$get_singleStr(private$combiStructure_index + 1)$get_seqFirstPos()
                  }
                },
                ## @description
                ## Private method: Get next singleStructureGenerator object in my_combiStructure object
                ## @return next singleStructureGenerator object if it exists, NULL if it does not exist
                get_nextStr = function(){
                  if (is.null(private$combiStructure_index)) return(NULL)
                  if (private$combiStructure_index == private$my_combiStructure$get_singleStr_number()) return(NULL)
                  if (private$combiStructure_index < private$my_combiStructure$get_singleStr_number()){
                    private$my_combiStructure$get_singleStr(private$combiStructure_index + 1)
                  }
                },
                ## @description
                ## Private method: Get previous singleStructureGenerator object in my_combiStructure object
                ## @return previous singleStructureGenerator object if it exists, NULL if it does not exist
                get_prevStr = function(){
                  if (is.null(private$combiStructure_index)) return(NULL)
                  if (private$combiStructure_index == 1) return(NULL)
                  if (private$combiStructure_index > 1){
                    private$my_combiStructure$get_singleStr(private$combiStructure_index - 1)
                  }
                },
                ## @description
                ## Private method: Update $neighbSt of a CpG position's neighbors within singleStructureGenerator object
                ##
                ## This fuction takes in the position index of a CpG site that
                ## changed $seq state and updates $neighbSt for the neighbors of the changed position
                ##
                ## @param position position index of the $seq change
                ##
                ## @return NULL
                ##
                update_intraStr_neighbSt = function(position){
                  if (!is.numeric(position) || length(position) != 1) {
                    stop("'position' must be one number")
                  }
                  if( position < 1 || position > length(private$seq)){
                    stop("'position' value must be within $seq length")
                  }

                  if (length(private$seq)== 1){ # cases with length 1
                    private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[1],private$seq[1]]
                  } else { # cases with length > 1
                    if (position == 1){
                      # position 1 has only one neighbor: position 2
                      if(length(private$seq) < position + 2){ # if position 2 has no right neighbor, it counts position 1 as both neighbors
                        private$neighbSt[position + 1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position]]
                      } else {
                        private$neighbSt[position + 1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position + 2]]
                      }
                    } else if (position == 2){
                      #1st counts the only neighbor as both neighbors
                      private$neighbSt[position-1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position]]
                      private$neighbSt[position+1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position+2]]
                    } else if (position == length(private$seq)){
                      # position n has only one neighbor
                      private$neighbSt[position-1] <<- private$mapNeighbSt_matrix[private$seq[position-2], private$seq[position]]
                    } else if (position == length(private$seq)-1){
                      #last counts the only neighbor as both neighbors
                      private$neighbSt[position+1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position]]
                      private$neighbSt[position-1] <<- private$mapNeighbSt_matrix[private$seq[position-2], private$seq[position]]
                    } else {
                      # intermediate positions
                      private$neighbSt[position-1] <<- private$mapNeighbSt_matrix[private$seq[position - 2], private$seq[position]]
                      private$neighbSt[position+1] <<- private$mapNeighbSt_matrix[private$seq[position], private$seq[position+2]]
                    }
                  }
                },
                ## @description
                ## Private method: Update $neighbSt of a CpG position's neighbors within and between singleStructureGenerator objects
                ##
                ## This fuction takes in the position index of a CpG site that
                ## changed $seq state and updates $neighbSt for the neighbors of the changed position
                ##
                ## @param position position index of the $seq change
                ##
                ## @return NULL
                ##
                update_neighbSt = function(position){
                  if (is.null(private$combiStructure_index)){ # case of isolated singleStructure instance
                    private$update_intraStr_neighbSt(position)
                  } else { # case of singleStructure instances within combiStructure

                    ## Update left neighbSt
                    if (position == 1){
                      if(!is.null(private$get_prevStr())){ # if singleStr is not first Str
                        if(!is.null(private$get_prevStr()$get_seq2ndButLastPos())){ # if there is 2 positions to the left
                          private$get_prevStr()$update_interStr_lastNeighbSt(private$get_prevStr()$get_seq2ndButLastPos(), private$seq[position])
                        } else {
                          private$get_prevStr()$update_interStr_lastNeighbSt(private$seq[position], private$seq[position])
                        }
                      }
                    } else if (position == 2){
                      if(!is.null(private$get_prevStr())){ # if singleStr is not first Str
                        private$neighbSt[position - 1] <- private$mapNeighbSt_matrix[private$get_leftStr_neighbSt(), private$seq[position]]
                      } else { # singleStr has no neighbouring Str to the left, first position takes as 2 neighbors second position
                        private$neighbSt[position - 1] <- private$mapNeighbSt_matrix[private$seq[position], private$seq[position]]
                      }
                    } else {
                      private$neighbSt[position - 1] <- private$mapNeighbSt_matrix[private$seq[position - 2], private$seq[position]]
                    }

                    ## Update right neighbSt
                    if (position == length(private$seq)){
                      if(!is.null(private$get_nextStr())){ # if singleStr is not last Str
                        if(!is.null(private$get_nextStr()$get_seq2ndPos())){ # if there is 2 positions to the right
                          private$get_nextStr()$update_interStr_firstNeighbSt(private$seq[position], private$get_nextStr()$get_seq2ndPos())
                        } else {
                          private$get_nextStr()$update_interStr_firstNeighbSt(private$seq[position], private$seq[position])
                        }
                      }
                    } else if (position == length(private$seq)-1){
                      if(!is.null(private$get_nextStr())){ # if singleStr is not last Str
                        private$neighbSt[position + 1] <- private$mapNeighbSt_matrix[private$seq[position], private$get_rightStr_neighbSt()]
                      } else { # singleStr has no neighbouring Str to the righ, last positions takes as 2 neighbors last - 1 position
                        private$neighbSt[position + 1] <- private$mapNeighbSt_matrix[private$seq[position], private$seq[position]]
                      }
                    } else {
                      private$neighbSt[position + 1] <- private$mapNeighbSt_matrix[private$seq[position], private$seq[position + 2]]
                    }
                  }
                },
                ## @field alpha_Ri Private attribute: Model parameter for gamma distribution shape to initialize the 3 $Ri_values
                alpha_Ri = 0.1,
                ## @field iota Private attribute: Model parameter for gamma distribution expected value to initialize the 3 $Ri_values
                iota = 0.3,
                ## @field Ri_values Private attribute: Vector containing 3 Ri rate values
                Ri_values = NULL,
                ## @description
                ## Private method: Initialize $Ri_values
                ##
                ## This fuction uses $iota and $alpha_Ri to compute the 3 Ri_values
                ##
                ## @return A list with the 3 Ri_values
                ##
                init_Ri_values = function(){
                  # Set the values that divide the gamma distribution in 3 equal probability categories
                  qGamma_oneThird <- stats::qgamma(1/3, shape = private$alpha_Ri, scale= private$iota / private$alpha_Ri)
                  qGamma_twoThird <- stats::qgamma(2/3, shape = private$alpha_Ri, scale= private$iota / private$alpha_Ri)

                  # Calculate the center of gravity of each gamma distribution category
                  Ri1 <- stats::integrate(function(x) x*stats::dgamma(x, shape =  alpha_Ri, scale =  iota / alpha_Ri) * 3, 0, qGamma_oneThird)$value
                  Ri2 <- stats::integrate(function(x) x*stats::dgamma(x, shape =  alpha_Ri, scale =  iota / alpha_Ri) * 3, qGamma_oneThird, qGamma_twoThird)$value
                  Ri3 <- stats::integrate(function(x) x*stats::dgamma(x, shape =  alpha_Ri, scale =  iota / alpha_Ri) * 3, qGamma_twoThird, Inf)$value
                  return(c(Ri1, Ri2, Ri3))
                },
                ## @field Rc_values Private attribute: Vector containing 2 Rc rate values
                Rc_values = NULL,
                ## @description
                ## Private method: Initialize $Rc_values
                ##
                ## This fuction uses $iota to compute the 2 Rc_values
                ##
                ## @return A list with the 2 Rc_values
                ##
                init_Rc_values = function(){
                  Rcl <- (1 - private$iota) / 2
                  Rcr <- (1 - private$iota) / 2
                  return(list(Rcl = Rcl, Rcr= Rcr))
                },
                ## @field Qi Private attribute: Rate matrix for the independent SSE process
                Qi = NULL,
                ## @description
                ## Private method: Set $Qi
                ##
                ## This fuction is used:
                ##
                ## 1. **During Initialization:** During the initialization of the object, this method is automatically called to set $Qi.
                ##
                ## 2. **In IWE_evol():** The method is called within the IWE_evol() to update $Qi after new $eqFreqs are sampled.
                ##
                ## @return NULL
                ##
                set_Qi = function(){
                  ## Extract methylation state frequencies
                  u <- private$eqFreqs[1]
                  p <- private$eqFreqs[2]
                  m <- private$eqFreqs[3]

                  ## Create the rate matrices for SSEind
                  ## Matrix without R
                  Qi0 <- matrix(c(-p-m, p, m, u, -u-m, m, u, p, -u-p), nrow = 3, byrow = TRUE)
                  ## Matrices with the 3 R values
                  Qi1 <- Qi0 * private$Ri_values[1]
                  Qi2 <- Qi0 * private$Ri_values[2]
                  Qi3 <- Qi0 * private$Ri_values[3]
                  ## Save rate matrices in a list, and validate them
                  private$Qi <<- list(Qi1, Qi2, Qi3)
                },
                ## @field Qc Private attribute: Rate matrix for the correlated SSE process
                Qc = NULL,
                ## @description
                ## Private method: Set $Qc
                ##
                ## This fuction is used to initialize the singleStructure object
                ##
                ## @return NULL
                ##
                set_Qc = function(){
                  ##### Matrices for SSEcor for neighbors state as coded in mapNeighbSt_matrix
                  Rcl <- private$Rc_values$Rcl
                  Rcr <- private$Rc_values$Rcr
                  #### #### uu = 1
                  Qc1 <- matrix(c(0, 0, 0, Rcl + Rcr, -Rcl - Rcr, 0, Rcl + Rcr, 0, -Rcl - Rcr), nrow = 3, byrow = TRUE)
                  #### #### up = 2
                  Qc2 <- matrix(c(-Rcr, Rcr, 0, Rcl , -Rcl, 0, Rcl, Rcr, -Rcl - Rcr), nrow = 3, byrow = TRUE)
                  #### #### um = 3
                  Qc3 <- matrix(c(-Rcr, 0, Rcr, Rcl , -Rcl - Rcr, Rcr, Rcl, 0, -Rcl), nrow = 3, byrow = TRUE)
                  #### #### pu = 4
                  Qc4 <- matrix(c(- Rcl, Rcl, 0, Rcr , -Rcr, 0, Rcr, Rcl, -Rcl - Rcr), nrow = 3, byrow = TRUE)
                  #### #### pp = 5
                  Qc5 <- matrix(c(-Rcl - Rcr, Rcl + Rcr, 0, 0, 0, 0, 0, Rcl + Rcr, -Rcl - Rcr), nrow = 3, byrow = TRUE)
                  #### #### pm = 6
                  Qc6 <- matrix(c(-Rcl - Rcr, Rcl, Rcr, 0, -Rcr, Rcr, 0, Rcl, -Rcl), nrow = 3, byrow = TRUE)
                  #### #### mu = 7
                  Qc7 <- matrix(c(-Rcl, 0, Rcl, Rcr , -Rcl - Rcr, Rcl, Rcr, 0, -Rcr), nrow = 3, byrow = TRUE)
                  #### #### mp = 8
                  Qc8 <- matrix(c(-Rcl - Rcr, Rcr, Rcl, 0, -Rcl, Rcl, 0, Rcr, -Rcr), nrow = 3, byrow = TRUE)
                  #### #### mm = 9
                  Qc9 <- matrix(c(-Rcl - Rcr, 0, Rcl + Rcr, 0, -Rcl - Rcr, Rcl + Rcr, 0, 0, 0), nrow = 3, byrow = TRUE)
                  private$Qc <<- list(Qc1, Qc2, Qc3, Qc4, Qc5, Qc6, Qc7, Qc8, Qc9)
                },
                ## @field Q Private attribute: Rate matrix for the complete SSE process
                Q = NULL,
                ## @description
                ## Private method: Set $Q
                ##
                ## This fuction is used:
                ##
                ## 1. **During Initialization:** During the initialization of the object, this method is automatically called to set $Q.
                ##
                ## 2. **In IWE_evol():** The method is called within the IWE_evol() to update $Q after new $Qi is updated with the new $eqFreqs.
                ##
                ## @return NULL
                ##
                set_Q = function() {
                  # Initialize a list to store the results
                  private$Q <<- list()

                  for (i in 1:length(private$Qi)) {
                    # Initialize a list to store the results for this Ri
                    Qi_results <- list()

                    for (j in 1:length(private$Qc)) {
                      # Sum the rate matrices element-wise
                      result_matrix <- private$Qi[[i]] + private$Qc[[j]]

                      # Add the result to the list for this Ri
                      Qi_results[[j]] <- result_matrix
                    }

                    # Add list of results for corresponding Ri to the result list
                    private$Q[[i]] <<- Qi_results

                  }
                },
                ## @field ratetree Private attribute: Cumulative Rate of Change Binary Tree
                ## This will contain a list of K vectors of length 1, 2, 4,.., N=2^(K-1), where
                ## 2^(K-1) >= length($seq) > 2^(K-2).
                ## The first length($seq) entries of ratetree[[K]][1:length($seq)] will be the rates
                ## of a change at that sequence positions, and the ratetree[[K]][length($seq):N]
                ## are all 0.
                ## ratetree[[1]] will sum(ratetree[[K]]) contain the total rate of changes in the sequence.
                ## ratetree[[2]] will contain the sum of the values in the first half of ratetree[[K]] and
                ## the sum of the values in the second half of ratetree[[K]].
                ## The four values in ratetree[[3]] will be the sums of the values in each quarter of the values
                ## in ratetree[[K]], and so on.
                ratetree = NULL,
                ## @description
                ## Private method: Update $ratetree after a CpG position's change of methylation state
                ##
                ## This fuction takes in the position index of a CpG site that
                ## changed $seq state and the new rate of change and updates $ratetree
                ##
                ## @param position position index of the $seq change
                ## @param rate new rate of change
                ##
                ## @return NULL
                ##
                update_ratetree = function(position, rate) {
                  ## check_ratetree("beginning of update_ratetree")
                  K <- length(private$ratetree)
                  private$ratetree[[K]][position] <<- rate
                  if (K > 1){
                    for(k in (K-1):1) {
                      position <- floor((position-1)/2)+1
                      private$ratetree[[k]][position] <<- sum(private$ratetree[[k+1]][c(2*position-1,
                                                                                        2*position)])
                    }
                  }
                  ## check_ratetree("end of update_ratetree")
                 },

                ## @description
                ## Private method: checkes whether $ratetree is consitent
                ##
                ## @param s a string that can be used to indicate where in the code the test failed
                ##
                check_ratetree = function(s) {
                  K <- length(private$ratetree)
                  if (K > 1){
                      for(k in (K-1):1) {
                          for(i in 1:length(private$ratetree[[k]])) {
                              if(abs(private$ratetree[[k]][i] - sum(private$ratetree[[k+1]][c(2*i-1, 2*i)])) >
                                 1e-8*max(private$ratetree[[k+1]][c(2*i-1, 2*i)])) {
                                  stop(paste(paste(private$ratetree,collapse="\n"),
                                             "multiregion ratetree inconsistent at", s, k, i))
                              }
                          }
                      }
                  }
                },

                ## @description
                ## Private method: Update $ratetree after a CpG position's change of methylation state
                ##
                ## This fuction takes in the position index of a CpG site that
                ## changed $seq state and the new rate of change and updates $ratetree
                ##
                ## @param testing default FALSE. TRUE for enabling testing output
                ##
                ## @return $seq index if testing FALSE and additional sampled value if testing TRUE
                ##
                choose_random_seqpos = function(testing=FALSE) {
                  ## return a random sequence position that is chosen with a probability that is
                  ## proportional to the change rates stored in private$ratetree
                  r  <- runif(n = 1, min = 0, max = private$ratetree[[1]][1])
                  s <- r
                  i <- 1
                  K <- length(private$ratetree)
                  if (K > 1){
                    for(k in 2:K) {
                      i <- 2*i-1
                      if(s > private$ratetree[[k]][i]) {
                        s <- s - private$ratetree[[k]][i]
                        i <- i+1
                      }
                    }
                  } else {
                    i = 1
                  }
                  if(testing){
                    return(c(i, r))
                  }
                  i
                },
                ## @description
                ## Private method: Sample a number of changes in time interval of length dt
                ##
                ## This fuction takes in the length of the time interval (dt) and
                ## samples a number of changing state events from a Poisson distribution
                ## with rate given by the total rate of change of the singleStructure object
                ## over dt time
                ##
                ## @param dt length of time interval
                ##
                ## @return sampled number of changes
                ##
                choose_number_of_changes = function(dt) {
                  ## choose a number of changes to happen in the next time interval of the short length dt
                  if (private$ratetree[[1]][1] < 0){
                    saveRDS(private$ratetree, "rateTreeError.rds")
                    rateMatrices <- list(Q = private$Q,
                                         Qi = private$Qi,
                                         Qc = private$Qc)
                    saveRDS(rateMatrices, "rateMatricesError.rds")
                    stop("negative total rate. Rate tree in file rateTreeError and rate matrices in rateMatricesError")
                  }
                  rpois(1, private$ratetree[[1]][1] * dt)
                }
              ),
              public = list(
                #' @description
                #' Public method: Initialization of $neighbSt
                #'
                #' This fuction initiates each CpG position $neighbSt as encoded in $mapNeighbSt_matrix
                #' Positions at the edge of the entire simulated sequence use
                #' their only neighbor as both neighbors.
                #'
                #' @return NULL
                init_neighbSt = function(){
                  if (is.null(private$my_combiStructure)){ # cases of singleStructure instances initiated outside combiStructure instance
                    if (length(private$seq)== 1){ # cases with length 1
                      private$neighbSt <<- private$mapNeighbSt_matrix[private$seq[1],private$seq[1]]
                    } else { # cases with length > 1
                      for (position in 1:length(private$seq)){
                        if (position == 1){ #1st counts the only neighbor as both neighbors
                          private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position + 1], private$seq[position + 1]]
                        } else if (position == length(private$seq)){ #last counts the only neighbor as both neighbors
                          private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position - 1]]
                        } else {
                          private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position + 1]]
                        }
                      }
                    }
                  } else { # cases of singleStructure instances initiated from combiStructure instance
                    if (length(private$seq)== 1){ # cases with length 1
                      if (private$combiStructure_index == 1){ # first singleStructure instance in combiStructure
                        if (is.null(private$get_nextStr())){ # fist singleStr is only singleStr
                          private$neighbSt <<- private$mapNeighbSt_matrix[private$seq[1],private$seq[1]]
                        } else { # first position next structure counts as both neighbors
                          private$neighbSt <<- private$mapNeighbSt_matrix[private$get_rightStr_neighbSt(), private$get_rightStr_neighbSt()]
                        }
                      } else if (private$combiStructure_index == private$my_combiStructure$get_singleStr_number()){ # last singleStructure instance in combiStructure, last position in previous structure counts as both neighbors
                        private$neighbSt <<- private$mapNeighbSt_matrix[private$get_leftStr_neighbSt(), private$get_leftStr_neighbSt()]
                      } else {
                        private$neighbSt <<- private$mapNeighbSt_matrix[private$get_leftStr_neighbSt(), private$get_rightStr_neighbSt()]
                      }
                    } else { # cases with length > 1
                      if (private$combiStructure_index == 1){ # first singleStructure instance in combiStructure
                        ## TODO: if length < position +2
                        for (position in 1:length(private$seq)){
                          if (position == 1){ #1st counts the only neighbor as both neighbors
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position + 1], private$seq[position + 1]]
                          } else if (position == length(private$seq)){ #last uses as right neighbor next singleStructure info
                            if(is.null(private$get_nextStr())){ # fist singleStr is only singleStr it counts previous position as both neighbors
                              private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position - 1]]
                            } else {
                              private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$get_rightStr_neighbSt()]
                            }
                          } else {
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position + 1]]
                          }
                        }
                      } else if (private$combiStructure_index == private$my_combiStructure$get_singleStr_number()){ # last singleStructure instance in combiStructure
                        for (position in 1:length(private$seq)){
                          if (position == 1){ #1st uses as left neighbor previous singleStructure info
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$get_leftStr_neighbSt(), private$seq[position + 1]]
                          } else if (position == length(private$seq)){ #last counts the only neighbor as both neighbors
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position - 1]]
                          } else {
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position + 1]]
                          }
                        }
                      } else { # intermediate singleStructure instances in combiStructure
                        for (position in 1:length(private$seq)){
                          if (position == 1){ #1st uses as left neighbor previous singleStructure info
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$get_leftStr_neighbSt(), private$seq[position + 1]]
                          } else if (position == length(private$seq)){ #last uses as right neighbor next singleStructure info
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$get_rightStr_neighbSt()]
                          } else {
                            private$neighbSt[position] <<- private$mapNeighbSt_matrix[private$seq[position - 1], private$seq[position + 1]]
                          }
                        }
                      }
                    }
                  }
                },
                #' @description
                #' Public method: Initialization of $ratetree
                #'
                #' This function initializes $ratetree
                #'
                #' @return NULL
                initialize_ratetree = function() {
                  K  <- ceiling(log(length(private$seq), 2)) + 1
                  for (k in 1:K) private$ratetree[[k]] <<- rep(0.0, 2^(k-1))
                  for (i in 1:length(private$seq)) {
                    private$ratetree[[K]][i] <<- abs(private$Q[[private$siteR[i]]][[private$neighbSt[i]]][private$seq[i],private$seq[i]])
                    ### AND NOTE THAT THE WAY THIS IS NOW CHOSEN MEANS THAT THE CHOSEN SITE SHOULD REALLY CHANGE (!) THE STATE
                  }
                  # Control for length of seq is 1, K = 1, k is 0
                  if(K > 1){
                    for (k in (K-1):1) {
                      for(i in 1:length(private$ratetree[[k]])) private$ratetree[[k]][i] <<- sum(private$ratetree[[k+1]][c(2*i-1, 2*i)])
                    }
                  }
                },
                #' @description
                #' Create a new singleStructureGenerator object.
                #'
                #' Note that this object is typically generated withing a combiStructureGenerator object.
                #'
                #' @param globalState Character. Structure's favored global state: "M" for methylated (island structures) / "U" for unmethylated (non-island structures).
                #' @param n Numerical Value. Number of CpG positions
                #' @param eqFreqs Default NULL. When given: numerical vector with structure's methylation state equilibrium frequencies (for unmethylated, partially methylated and methylated)
                #' @param combiStr Default NULL. When initiated from combiStructureGenerator: object of class combiStructureGenerator containing it
                #' @param combiStr_index Default NULL. When initiated from combiStructureGenerator: index in Object of class combiStructureGenerator
                #' @param params Default NULL. When given: data frame containing model parameters
                #' @param testing Default FALSE. TRUE for testing output
                #' @return A new `singleStructureGenerator` object.
                initialize = function(globalState, n, eqFreqs = NULL, combiStr = NULL, combiStr_index = NULL,  params = NULL, testing = FALSE) {
                  if (!is.character(globalState)) {
                    stop("globalState must be a character")
                  }
                  # Check if globalState is valid
                  if (!(globalState %in% c("M", "U"))) {
                    stop("globalState must be 'M' or 'U'")
                  }
                  if (!is.numeric(n)) {
                    stop("n must be a numeric value")
                  }
                  if (!length(n)==1){
                    stop("n must be of length 1")
                  }
                  private$globalState <- globalState
                  if(is.null(eqFreqs)){
                    private$eqFreqs <- private$sample_eqFreqs()
                  } else {
                    if(!is.numeric(eqFreqs) || !length(eqFreqs) == 3 || !sum(eqFreqs)==1){
                      stop("if 'eqFreqs' is not NULL, provide a numeric vector of 3 frequencies ")
                    }
                    private$eqFreqs <- eqFreqs
                  }
                  if(!is.null(params)){
                    private$alpha_pI <- params$alpha_pI
                    private$beta_pI <- params$beta_pI
                    private$alpha_mI <- params$alpha_mI
                    private$beta_mI <- params$beta_mI
                    private$alpha_pNI <- params$beta_pNI
                    private$beta_pNI <- params$beta_pNI
                    private$alpha_mNI <- params$alpha_mNI
                    private$beta_mNI <- params$beta_mNI
                    private$alpha_Ri <- params$alpha_Ri
                    private$iota <- params$iota
                  }
                  private$seq <- sample(1L:3L, size = n, prob = private$eqFreqs, replace = TRUE)
                  if(testing){ # when testing neighbSt initiate instance with n=13
                    private$seq <- c(1, 1, 2, 3, 1, 1, 1, 3, 2, 2, 3, 2, 3)
                  }
                  private$siteR <- sample(1L:3L, size = n, replace = TRUE)
                  private$my_combiStructure <- combiStr
                  private$combiStructure_index <- combiStr_index
                  private$Ri_values <- private$init_Ri_values()

                  private$Rc_values <- private$init_Rc_values()
                  private$set_Qi()
                  private$set_Qc()
                  private$set_Q()
                  #debug(self$init_neighbSt)
                  if(is.null(private$my_combiStructure)){
                    self$init_neighbSt()
                    #debug(self$initialize_ratetree)
                    self$initialize_ratetree()
                  }
                },
                #' @description
                #' Public method: Get object's methylation state sequence
                #'
                #' Encoded with 1 for unmethylated, 2 for partially methylated and 3 for methylated
                #'
                #' @return vector with equilibrium frequencies of unmethylated, partially methylated and methylated
                get_seq = function() private$seq,
                #' @description
                #' Public method: Get first sequence position methylation state
                #'
                #' @return numerical encoding of first position's methylation state
                get_seqFirstPos = function () private$seq[1],
                #' @description
                #' Public method: Get second sequence position methylation state
                #'
                #' @return numerical encoding of second position's methylation state. NULL if position does not exist
                get_seq2ndPos = function () {
                  if (length(private$seq)==1){
                    if (is.null(private$get_nextStr())) NULL
                    else private$get_nextStr()$get_seqFirstPos()
                  } else {
                    private$seq[2]
                  }
                },
                #' @description
                #' Public method: Get first sequence position methylation state
                #'
                #' @return numerical encoding of first position's methylation state
                get_seqLastPos = function() private$seq[length(private$seq)],
                #' @description
                #' Public method: Get second but last sequence position methylation state
                #'
                #' @return numerical encoding of second but last position's methylation state. NULL if position does not exist
                get_seq2ndButLastPos = function () {
                  if (length(private$seq)==1){
                    if (is.null(private$get_prevStr())) NULL
                    else private$get_prevStr()$get_seqLastPos()
                  } else {
                    private$seq[length(private$seq)-1]
                  }
                },
                #' @description
                #' Public method: Get index in object of class combiStructureGenerator
                #'
                #' @return index in object of class combiStructureGenerator
                get_combiStructure_index = function() private$combiStructure_index,
                #' @description
                #' Public method: Update neighbSt of next singleStructureGenerator object within combiStructureGenerator object
                #'
                #' This function is used when the last $seq position of a singleStructureGenerator object
                #' changes methylation state to update the neighbSt position
                #'
                #' @param leftNeighbSt $seq state of left neighbor (left neighbor is in previous singleStructureGenerator object)
                #' @param rightNeighbSt $seq state of right neighbor
                #'
                #' @return NULL
                update_interStr_firstNeighbSt = function(leftNeighbSt, rightNeighbSt) {
                  if(!is.null(rightNeighbSt)){
                    private$neighbSt[1] <<- private$mapNeighbSt_matrix[leftNeighbSt, rightNeighbSt]
                  } else { # if there is no right neighbour
                    private$neighbSt[1] <<- private$mapNeighbSt_matrix[leftNeighbSt, leftNeighbSt]
                  }
                },
                #' @description
                #' Public method: Update neighbSt of previous singleStructureGenerator object within combiStructureGenerator object
                #'
                #' @param leftNeighbSt $seq state of right neighbor (left neighbor is in next singleStructureGenerator object)
                #' @param rightNeighbSt $seq state of right neighbor
                #'
                #' @return NULL
                update_interStr_lastNeighbSt = function(leftNeighbSt, rightNeighbSt){
                  if(!is.null(leftNeighbSt)){
                    private$neighbSt[length(private$neighbSt)] <<- private$mapNeighbSt_matrix[leftNeighbSt, rightNeighbSt]
                  } else { # if there is no left neighbour
                    private$neighbSt[length(private$neighbSt)] <<- private$mapNeighbSt_matrix[rightNeighbSt, rightNeighbSt]
                  }
                },
                #' @description
                #' Public method: Get object's equilibrium Frequencies
                #'
                #' @return vector with equilibrium frequencies of unmethylated, partially methylated and methylated
                get_eqFreqs = function() private$eqFreqs,

                #' @description
                #' Public method. Simulate how CpG dinucleotide methylation state changes due to the SSE process
                #' along a time step of length dt
                #'
                #' @param dt time step length.
                #' @param testing logical value for testing purposes. Default FALSE.
                #'
                #' @return default NULL. If testing TRUE it returns a list with the number of events sampled and a
                #' dataframe with the position(s) affected, new state and old methylation state.
                #'
                SSE_evol = function(dt, testing = FALSE) {
                    if (testing){
                        event_number <- 0
                        SSE_evolInfo <- data.frame(
                            position = integer(),
                            old_St = integer(),
                            new_St = integer()
                        )
                    }
                    ## get a number of changes to happen in the next time interval of the short length dt
                    M <- private$choose_number_of_changes(dt)
                    if (M>0){
                        for(m in 1:M) {
                            i <- private$choose_random_seqpos()
                            if (testing){
                                event_number <- M
                                position <- i
                                old_St <- private$seq[i]
                            }
                                        # assign new sequence position state with probability given by the relative rates of changing to each of the 2 other states
                            # Define a global variable to store error information
                            error_info_global <- NULL
                            
                            # Try-catch block with error handling
                            error_info <- tryCatch({
                              # Attempt to run sample()
                              private$seq[i] <<- sample(1:3, size = 1, prob = sapply(Q[[private$siteR[i]]][[private$neighbSt[i]]][private$seq[i], ], max, 0))
                              NULL # If no error, return NULL
                            }, error = function(e) {
                              # Capture error information in a list
                              singleStrcloned <- self$clone()
                              list(
                                message = conditionMessage(e),
                                ratetree = private$ratetree[[1]][1],
                                M = M,
                                probabilities = sapply(Q[[private$siteR[i]]][[private$neighbSt[i]]][private$seq[i], ], max, 0),
                                Q = Q[[private$siteR[i]]][[private$neighbSt[i]]][private$seq[i], ],
                                singleStructure_index = private$combiStructure_index,
                                CpG_position = i,
                                singleStrcloned = singleStrcloned
                              )
                            })
                            
                            # If error_info is not NULL, it means an error occurred
                            if (!is.null(error_info)) {
                              print("Error occurred:")
                              print(error_info$message)
                              print(paste("ratetree", error_info$ratetree))
                              print(paste("M", error_info$M))
                              print("Probabilities")
                              print(error_info$probabilities)
                              print("Q")
                              print(error_info$Q)
                              print(paste("my_singleStructure_index", error_info$singleStructure_index))
                              print(paste("CpG position", error_info$CpG_position))
                              # Print singleStrcloned information if needed
                              print("singleStrcloned")
                              print(error_info$singleStrcloned)
                              
                              # Save the error information to a global variable
                              error_info_global <<- error_info
                              
                              # Stop execution
                              stop("Execution stopped due to an error.")
                            }
                            
                            
                            
                            
                            #private$seq[i] <<- sample(1:3, size=1, prob=sapply(Q[[private$siteR[i]]][[private$neighbSt[i]]][private$seq[i],], max, 0))
                            if (testing){
                                new_St <- private$seq[i]
                                SSE_evolInfo <- rbind(SSE_evolInfo, data.frame(position, old_St, new_St))
                            }
                            private$update_neighbSt(i)
                            for(j in max(i-1, 1):min(i+1, length(private$seq))) {
                                private$update_ratetree(j, abs(private$Q[[private$siteR[j]]][[private$neighbSt[j]]][private$seq[j],private$seq[j]]))
                            }
                        }
                    }
                    if (testing){
                        return(list(event_number = event_number,
                                    SSE_evolInfo = SSE_evolInfo))
                    }
                },

                #' @description
                #' Public Method. Simulate IWE Events
                #'
                #' Simulates how CpG Islands' methylation state frequencies change and simultaneous sites change methylation state
                #' along a branch of length t according to the SSE-IWE model.
                #'
                #' @param testing logical value for testing purposes. Default FALSE.
                #'
                #' @return If testing = TRUE it returns a list.
                #' If there was a change in the equilibrium frequencies the list contains the following 7 elements, if not it contains the first 3 elements:
                #' \describe{
                #'       \item{\code{eqFreqsChange}}{logical indicating if there was a change in the equilibrium frequencies.}
                #'       \item{\code{old_eqFreqs}}{Original equilibrium frequencies before the IWE event.}
                #'       \item{\code{new_eqFreqs}}{New equilibrium frequencies after the IWE event.}
                #'       \item{\code{old_obsFreqs}}{Original observed frequencies before the IWE event.}
                #'       \item{\code{new_obsFreqs}}{New observed frequencies after the IWE event.}
                #'       \item{\code{IWE_case}}{Description of the IWE event case.}
                #'       \item{\code{Mk}}{Transition matrix used for the IWE event.}
                #' }
                #'
                #'
                #' @details The function checks if the methylation equilibrium frequencies (\code{eqFreqs}) and sequence observed
                #' frequencies (\code{obsFreqs}) change after the IWE event. If there is a change in either
                #' frequencies, the corresponding change flags(\code{eqFreqsChange}
                #' in the \code{infoIWE} list will be set to \code{TRUE}.
                #'
                IWE_evol = function(testing = FALSE) {
                    ## Extract previous state equilibrium frequencies
                    u <- private$eqFreqs[1]
                    p <- private$eqFreqs[2]
                    m <- private$eqFreqs[3]
                    old_eqFreqs <- c(u,p,m)

                    if (testing){
                                        # compute previous observed methylation frequencies
                        old_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
                                        # Save previous rates
                        old_Q <- private$Q
                                        # Initiate changedPos with NULL
                        changedPos <- NULL
                    }

                                        # Sample new methylation frequencies with the structure's global strategie
                    new_eqFreqs <- private$sample_eqFreqs()

                                        # Check if new_eqFreqs are equal to old ones
                    if (identical(old_eqFreqs, new_eqFreqs)) {
                        eqFreqsChange = F

                    } else{
                        eqFreqsChange = T

                                        # Set the IWE transition matrix according to current case
                                        # Check Case 1: 1 new frequency value bigger and 2 smaller
                        if (new_eqFreqs[1] > u & new_eqFreqs[2] <= p & new_eqFreqs[3] <= m) {
                            IWE_case <- "Case 1. u bigger"
                            Mk <- matrix(c(1, 0, 0,
                            (p-new_eqFreqs[2])/p, new_eqFreqs[2]/p, 0,
                            (m-new_eqFreqs[3])/m, 0, new_eqFreqs[3]/m),
                            nrow = 3, byrow = TRUE)

                        }
                        if (new_eqFreqs[2] > p & new_eqFreqs[1] <= u & new_eqFreqs[3] <= m) {
                            IWE_case <- "Case 1. p bigger"
                            Mk <- matrix(c(new_eqFreqs[1]/u, (u-new_eqFreqs[1])/u, 0,
                                           0, 1, 0,
                                           0, (m-new_eqFreqs[3])/m, new_eqFreqs[3]/m),
                                         nrow = 3, byrow = TRUE)
                        }
                        if (new_eqFreqs[3] > m & new_eqFreqs[2] <= p & new_eqFreqs[1] <= u) {
                            IWE_case <- "Case 1. m bigger"
                            Mk <- matrix(c(new_eqFreqs[1]/u, 0, (u-new_eqFreqs[1])/u,
                                           0, new_eqFreqs[2]/p, (p-new_eqFreqs[2])/p,
                                           0, 0, 1),
                                         nrow = 3, byrow = TRUE)

                        }
                                        # Check Case 2: 1 new frequency value smaller
                        if (new_eqFreqs[1] < u & new_eqFreqs[2] >= p & new_eqFreqs[3] >= m) {
                            IWE_case <- "Case 2. u smaller"
                            Mk <- matrix(c(new_eqFreqs[1]/u, (new_eqFreqs[2]-p)/u, (new_eqFreqs[3]-m)/u,
                                           0, 1, 0,
                                           0, 0, 1),
                                         nrow = 3, byrow = TRUE)
                        }
                        if (new_eqFreqs[2] < p & new_eqFreqs[1] >= u & new_eqFreqs[3] >= m) {
                            IWE_case <- "Case 2. p smaller"
                            Mk <- matrix(c(1, 0, 0,
                            (new_eqFreqs[1]-u)/p, new_eqFreqs[2]/p, (new_eqFreqs[3]-m)/p,
                            0, 0, 1),
                            nrow = 3, byrow = TRUE)
                        }
                        if (new_eqFreqs[3] < m & new_eqFreqs[2] >= p & new_eqFreqs[1] >= u) {
                            IWE_case <- "Case 2. m smaller"
                            Mk <- matrix(c(1, 0, 0,
                                           0, 1, 0,
                                           (new_eqFreqs[1]-u)/m, (new_eqFreqs[2]-p)/m, new_eqFreqs[3]/m),
                                         nrow = 3, byrow = TRUE)
                        }

                                        # Change data equilibrium frequencies
                        private$eqFreqs <<- new_eqFreqs
                                        # Update Qi and Q
                        private$set_Qi()
                        private$set_Q()

                        if(testing){

                                        # Validate updated rate matrices
                            validationStates <- listRateMatrix_validation(private$Qi, "Qi matrices IWE")
                            listMatrices_validationResults(validationStates)

                            validationStates <- listRateMatrix_validation(private$Q[[1]], "Q matrices IWE R1")
                            listMatrices_validationResults(validationStates)

                            validationStates <- listRateMatrix_validation(private$Q[[2]], "Q matrices IWE R2")
                            listMatrices_validationResults(validationStates)

                            validationStates <- listRateMatrix_validation(private$Q[[3]], "Q matrices IWE R3")
                            listMatrices_validationResults(validationStates)

                                        # Validate transition matrix
                            validationStates <- listTransitionMatrix_validation(list(Mk), "Mk")
                            listMatrices_validationResults(validationStates)

                                        # Validate methylation equilibrium tripple (extra control at IWE)
                            validationStates <- listFreqVector_validation(list(old_eqFreqs, new_eqFreqs), "IWE_evol control old_eqFreqs and new_eqFreqs")
                            listFreqVector_validationResults(validationStates)

                                        # Validate Markov Chain State Transition Property
                            validationStates <- transPropMC_validation(old_eqFreqs, Mk, new_eqFreqs, listName = "transPropMC IWE")
                            transPropMC_validationResults(validationStates)
                        }
                        
                        if(!exists("Mk")){
                          print("new_eqFreqs")
                          print(new_eqFreqs)
                          print("old_eqFreqs")
                          print(old_eqFreqs)
                        }
                        #print("in IWE")

                                        # Sample $seq accordint to transition probablities
                        newseq <- rep(0, length(private$seq))
                        for(i in 1:length(newseq)) {
                            newseq[i] <- sample(1:3, size=1, prob=as.vector(Mk[private$seq[i],]))
                        }

                                        # Update $seq and $neighbSt
                        if(any(private$seq != newseq)){
                            changedPos <- which(private$seq !=newseq)
                            private$seq <<- newseq
                            for (i in changedPos){
                                private$update_neighbSt(i)
                                        # Update $ratetree
                                #for(j in max(i-1, 1):min(i+1, length(private$seq))) {
                                #    private$update_ratetree(j, abs(private$Q[[private$siteR[j]]][[private$neighbSt[j]]][private$seq[j],private$seq[j]]))
                                #}
                            }
                        }
                        self$initialize_ratetree()

                                        #Compute the new_obsFreqs after the IWE event
                        new_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
                    }

                    if (testing){
                        if(eqFreqsChange){
                            return(list(combiStructure_index = private$combiStructure_index,
                                        eqFreqsChange = eqFreqsChange,
                                        old_eqFreqs = old_eqFreqs,
                                        new_eqFreqs = new_eqFreqs,
                                        old_obsFreqs = old_obsFreqs,
                                        new_obsFreqs = new_obsFreqs,
                                        IWE_case = IWE_case,
                                        Mk = Mk,
                                        old_Q = old_Q,
                                        new_Q = private$Q,
                                        changedPos = changedPos))
                        } else {
                            return(list(eqFreqsChange = eqFreqsChange,
                                        old_eqFreqs = old_eqFreqs,
                                        new_eqFreqs = new_eqFreqs))
                        }
                    }
                },
                                        # Island equilibrium frequencies sampling strategia
                #' @description
                #' Public Method.
                #' @return Model parameter alpha_pI for sampling island equilibrium frequencies
                get_alpha_pI = function() private$alpha_pI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling island equilibrium frequencies
                get_beta_pI = function() private$beta_pI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling island equilibrium frequencies
                get_alpha_mI = function() private$alpha_mI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling island equilibrium frequencies
                get_beta_mI = function() private$beta_mI,

                                        #sampling strategy
                #' @description
                #' Public Method.
                #' @return Model parameter for sampling non-island equilibrium frequencies
                get_alpha_pNI = function() private$alpha_pNI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling non-island equilibrium frequencies
                get_beta_pNI = function() private$beta_pNI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling non-island equilibrium frequencies
                get_alpha_mNI = function() private$alpha_mNI,

                #' @description
                #' Public Method.
                #' @return Model parameter for sampling non-island equilibrium frequencies
                get_beta_mNI = function() private$beta_mNI,

                                        # SSE
                #' @description
                #' Public Method.
                #' @return Model parameter for gamma distribution shape to initialize the 3 $Ri_values
                get_alpha_Ri = function() private$alpha_Ri,

                #' @description
                #' Public Method.
                #' @return Model parameter for gamma distribution expected value to initialize the 3 $Ri_values
                get_iota = function() private$iota,

                #' @description
                #' Public Method.
                #' @return The 3 $Ri_values
                get_Ri_values = function() private$Ri_values
              )
              )

#' @title combiStructureGenerator
#' @importFrom R6 R6Class
#'
#' @description
#' an R6 class representing several genomic structures.
#' Each genomic structure contained is an object of class singleStructureGenerator.
#' Note that default clone(deep=TRUE) fails to clone singleStructureGenerator objects contained, use method $copy() instead.
#'
combiStructureGenerator <-
  R6::R6Class("combiStructureGenerator",
              private = list(
                ## @field singleStr Private attribute: List containing objects of class singleStructureGenerator
                singleStr = NULL,
                ## @field globalState Private attribute: Vector with each singleStructureGenerator object favored global state
                singleStr_globalState = NULL,
                ## @field mu Private attribute: Model parameter for the rate of the IWE evolutionary process (per island and branch length).
                mu = 0.1,
                ## @field IWE_rate Private attribute: Total rate of the IWE evolutionary process for the combiStructureObject.
                IWE_rate = NULL,
                ## @description
                ## Private Method: Compute the total rate of IWE evolutionary process.
                ## The total rate of change is given by the rate per CpG island times the number of CpG islands
                ##
                ## @return NULL.
                set_IWE_rate = function(){
                  if (self$get_island_number() != 0){
                    private$IWE_rate <- self$get_island_number() * private$mu
                  } else {
                    private$IWE_rate <- 0
                  }
                },

                ## @description
                ## Private method: Simulate SSE evolution at each short time step in each of the
                ## singleStructureGenerator objects in $singleStr. It samples each object in $singleStr
                ## in random order and calls singleStructureGenerator$SSE_evol()
                ##
                ## @param dt length of time step for SSE evolution.
                ## @param testing default FALSE. True for testing output
                ##
                ## @return default NULL. If testing TRUE it returns testing information
                SSE_evol = function(dt, testing = FALSE){
                    evol_order <- sample(1:self$get_singleStr_number(), self$get_singleStr_number(), replace = FALSE)
                    if (!testing){
                        for (i in evol_order) private$singleStr[[i]]$SSE_evol(dt, testing)
                    } else {
                        testing_info <- vector("list", self$get_singleStr_number())
                        for (i in evol_order) testing_info[[i]] <- private$singleStr[[i]]$SSE_evol(dt, testing)
                    }
                    if (testing){
                        return(list(evol_order = evol_order,
                                    testing_info = testing_info))
                    }
                },

                ## @description
                ## Simulate SSE evolutionary process within a given interval length by time steps of length dt.
                ##
                ## @param testing Default FALSE. Logical for testing purposes.
                ## @param dt Length of time steps for SSE_evol.
                ## @param interval_length Length of the total time interval.
                ##
                ## @return Default NULL. If testing = TRUE returns information for testing purposes.
                ##
                ## @details The function simulates SSE events within a specified interval length.
                ##
                interval_evol = function(interval_length, dt, testing = FALSE) {
                    if(testing){
                        interval_evolInfo <- c()
                    }
                                        # Compute number of dt for discretized time steps
                    number_dt <- floor(interval_length / dt)
                    remainder <- interval_length %% dt
                    if(interval_length < dt){
                        if(testing){
                            interval_evolInfo <- paste("SSE evolving in shorter interval_length than dt. Interval length:", interval_length)
                        }
                        private$SSE_evol(dt = interval_length)
                    } else{
                        if(testing){
                            interval_evolInfo <- paste("Number of SSE:", number_dt, ". Interval length:", dt)
                        }
                        for(i in 1:number_dt){
                            private$SSE_evol(dt = dt)
                        }
                        if (dt * number_dt == interval_length){
                            if(testing){
                                interval_evolInfo <- c(interval_evolInfo, "No remainder")
                            }
                        } else{
                            if(testing){
                                interval_evolInfo <- c(interval_evolInfo, paste("Last SSE with remainder dt of length", remainder))
                            }
                            private$SSE_evol(dt = remainder)
                        }
                    }
                    if(testing){interval_evolInfo}
                },

                ## @field IWE_events Private attribute: Information of the IWE events sampled in a tree branch
                IWE_events = NULL,
                ## @field name Private attribute: If evolutionary process (simulated from class treeMultiRegionSimulator) ends in a tree leaf, the name of the leaf
                name = NULL,
                ## @field own_index Private attribute: Own branch index in the tree along which the evolutionary process is simulated (from class treeMultiRegionSimulator)
                own_index = NULL,
                ## @field parent_index Private attribute: Parent branch index in the tree along which the evolutionary process is simulated (from class treeMultiRegionSimulator)
                parent_index = NULL,
                ## @field parent_index Private attribute: Offspring branch index in the tree along which the evolutionary process is simulated (from class treeMultiRegionSimulator)
                offspring_index = NULL
              ),
              public = list(
                  #' @description
                  #' Create a new combiStructureGenerator object.
                  #'
                  #' Note that this object can be generated within a treeMultiRegionSimulator object.
                  #'
                  #' @param infoStr A data frame containing columns 'n' for the number of sites, and 'globalState' for the favoured global methylation state.
                  #' If initial equilibrium frequencies are given the dataframe must contain 3 additional columns: 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq'
                  #' @param params Default NULL. When given: data frame containing model parameters.
                  #' @param testing Default FALSE. TRUE for testing output.
                  #'
                  #' @return A new `combiStructureGenerator` object.
                  initialize = function (infoStr, params = NULL, testing = FALSE){
                      private$singleStr <- list()
                      private$singleStr_globalState <- c()
                      if (testing){ # data with seqlength 13
                          for (i in 1:nrow(infoStr)) {
                              u_length <- infoStr[i, "n"]
                              private$singleStr_globalState[i] <- infoStr[i, "globalState"]
                              private$singleStr[[i]] <- singleStructureGenerator$new(private$singleStr_globalState[i],
                                                                                     u_length,
                                                                                     combiStr = self, combiStr_index = i,
                                                                                     testing = TRUE)
                          }
                      } else {
                          for (i in 1:nrow(infoStr)) {
                              u_length <- infoStr[i, "n"]
                              private$singleStr_globalState[i] <- infoStr[i, "globalState"]
                              if(all(c("u_eqFreq", "p_eqFreq", "m_eqFreq") %in% colnames(infoStr))){
                                  eqFreqs <- c(infoStr$u_eqFreq[i], infoStr$p_eqFreq[i], infoStr$m_eqFreq[i])
                              } else{
                                  eqFreqs <- NULL
                              }
                              private$singleStr[[i]] <- singleStructureGenerator$new(globalState = private$singleStr_globalState[i],
                                                                                     n = u_length,
                                                                                     eqFreqs = eqFreqs,
                                                                                     combiStr = self ,combiStr_index = i,
                                                                                     params = params)

                          }
                      }
                      for (i in 1:length(private$singleStr)){
                          private$singleStr[[i]]$init_neighbSt()
                          private$singleStr[[i]]$initialize_ratetree()
                      }
                      if(!is.null(params)){
                          private$mu <- params$mu
                      }
                      private$set_IWE_rate()

                  },
                  #' @description
                  #' Public method: Get one singleStructureGenerator object in $singleStr
                  #'
                  #' @param i index of the singleStructureGenerator object in $singleStr
                  #'
                  #' @return the singleStructureGenerator object in $singleStr with index i
                  get_singleStr = function(i) private$singleStr[[i]],


                  #' @description
                  #' Public method: Get number of singleStructureGenerator objects in $singleStr
                  #'
                  #' @return number of singleStructureGenerator object contained in $singleStr
                  get_singleStr_number = function () length(private$singleStr),


                  #' @description
                  #' Public method: Get number of singleStructureGenerator objects in $singleStr with $globalState "U" (CpG islands)
                  #'
                  #' @return number of singleStructureGenerator in $singleStr objects with $globalState "U" (CpG islands)
                  get_island_number = function (){
                      indexes <- which(private$singleStr_globalState == "U")
                  if (length(indexes) == 0){
                    return(0)
                  } else {
                    return(length(indexes))
                  }
                },

                #' @description
                #' Public method: Get index of singleStructureGenerator objects in $singleStr with $globalState "U" (CpG islands)
                #'
                #' @return index of singleStructureGenerator objects in $singleStr with $globalState "U" (CpG islands)
                get_island_index = function() which(private$singleStr_globalState == "U"),


                #' @description
                #' Public method: Set information of the IWE events sampled in a tree branch
                #'
                #' @param a value to which IWE_events should be set
                #'
                #' @return NULL
                set_IWE_events = function(a) private$IWE_events <- a,


                #' @description
                #' Public method: Get information of the IWE events sampled in a tree branch
                #'
                #' @return information of the IWE events sampled in a tree branch
                get_IWE_events = function() private$IWE_events,


                #' @description
                #' Public method: Set the name of the leaf if evolutionary process
                #' (simulated from class treeMultiRegionSimulator) ends in a tree leaf.
                #'
                #' @param a value to which name should be set
                #'
                #' @return NULL
                set_name = function(a) {private$name <- a},


                #' @description
                #' Public method: Get the name of the leaf if evolutionary process
                #' (simulated from class treeMultiRegionSimulator) ended in a tree leaf.
                #'
                #' @return Name of the leaf if evolutionary process
                #' (simulated from class treeMultiRegionSimulator) ended in a tree leaf.
                #' For iner tree nodes return NULL
                get_name = function() {private$name},


                #' @description
                #' Public method: Set own branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @return NULL
                get_own_index =  function() {private$own_index},


                #' @description
                #' Public method: Get own branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @param i index of focal object
                #'
                #' @return Own branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                set_own_index =  function(i) {private$own_index <- i},


                #' @description
                #' Public method: Get parent branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @return Parent branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                get_parent_index = function() {private$parent_index},


                #' @description
                #' Public method: Set parent branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @param i set parent_index to this value
                #'
                #' @return NULL
                set_parent_index = function(i) {private$parent_index <- i},


                #' @description
                #' Public method: Get offspring branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @return Offspring branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                get_offspring_index = function() {private$offspring_index},


                #' @description
                #' Public method: Set offspring branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @param i set offspring_index to this value
                #'
                #' @return NULL
                set_offspring_index = function(i) {private$offspring_index <- i},


                #' @description
                #' Public method: Add offspring branch index in the tree
                #' along which the evolutionary process is simulated
                #' (from class treeMultiRegionSimulator).
                #'
                #' @param i index to be added
                #'
                #' @return NULL
                add_offspring_index = function(i) {
                  if(is.null(private$offspring_index)) {
                    private$offspring_index <- i
                  } else {
                    private$offspring_index <- c(private$offspring_index, i)
                  }
                },

                #' @description
                #' Public method.
                #'
                #' @return Model parameter for the rate of the IWE evolutionary process (per island and branch length).
                get_mu = function() private$mu,

                #' @description
                #' Public method: Clone each singleStructureGenerator object in $singleStr
                #'
                #' @param singStrList object to be cloned
                #'
                #' @return NULL
                set_singleStr = function(singStrList){
                    private$singleStr <- lapply(singStrList, function(singleStr) {
                        singleStr$clone()
                    })
                },

                #' @description
                #' Public method: Clone combiStructureGenerator object and all singleStructureGenerator objects in it
                #'
                #' @return cloned combiStructureGenerator object
                copy = function(){
                    new_obj <- self$clone()
                    new_obj$set_singleStr(private$singleStr)
                    return(new_obj)
                },

                #' @description
                #' Simulate CpG dinucleotide methylation state evolution along a tree branch.
                #' The function samples the IWE events on the tree branch and simulates the
                #' evolution through the SSE and IWE processes.
                #'
                #' @param dt Length of SSE time steps.
                #' @param testing Default FALSE. TRUE for testing purposes.
                #' @param branch_length Length of the branch.
                #'
                #' @return Default NULL. If testing = TRUE it returns information for testing purposes.
                #'
                #' @details
                #' It handles both cases where IWE events are sampled or not sampled within the branch.
                #'
                branch_evol = function(branch_length, dt, testing = FALSE) {
                    branchEvolInfo <- list()
                                        # Initialize IWE_times as an empty vector
                                        #branchEvolInfo$IWE_times <- c()

                    if(private$IWE_rate == 0){
                        branchEvolInfo$IWE_event <- FALSE
                        self$set_IWE_events("Simulation without IWE events.")
                        SSE_intervals <- branch_length
                                        # Sequence evolution along the branch
                        private$interval_evol(interval_length = branch_length, dt = dt)
                    } else {
                                        # Draw IWE times
                        IWE_t <- stats::rexp(1, private$IWE_rate)
                        if (IWE_t >= branch_length) { # If there is no IWE sampled in the tree branch
                            branchEvolInfo$IWE_event <- FALSE
                                        # Save IWE_events
                            self$set_IWE_events(FALSE)
                            branchEvolInfo$IWE_times <- IWE_t
                            SSE_intervals <- branch_length
                                        # Sequence evolution along the branch
                            private$interval_evol(interval_length = branch_length, dt = dt)
                        } else{ # If there is at least 1 IWE sampled in the tree branch
                            branchEvolInfo$IWE_event <- TRUE
                            while (IWE_t < branch_length) {
                                branchEvolInfo$IWE_times <- c(branchEvolInfo$IWE_times, IWE_t)
                                IWE_t <- IWE_t + stats::rexp(1, private$IWE_rate)
                            }
                                        # Sample islands to apply IWE_events
                            if(length(self$get_island_index()) == 1){
                                branchEvolInfo$islands <- rep(self$get_island_index(), length(branchEvolInfo$IWE_times))
                            } else {
                                branchEvolInfo$islands <- sample(self$get_island_index(), length(branchEvolInfo$IWE_times), replace = TRUE)
                            }
                                        # Save islands and IWE times
                            self$set_IWE_events(list(islands = branchEvolInfo$islands,
                                                     times = branchEvolInfo$IWE_times))
                                        # Compute the branch intervals for SSE evolution
                            if (length(branchEvolInfo$IWE_times)==1){
                                SSE_intervals <- c(branchEvolInfo$IWE_times, branch_length-branchEvolInfo$IWE_times)
                            } else if (length(branchEvolInfo$IWE_times) > 1){
                                SSE_intervals <- c(branchEvolInfo$IWE_times[1], diff(branchEvolInfo$IWE_times), branch_length - branchEvolInfo$IWE_times[length(branchEvolInfo$IWE_times)])
                            }
                            branchEvolInfo$SSE_intervals <- SSE_intervals
                            branchEvolInfo$infoIWE <- list()
                                        # Sequence evolution along the branch intervals followed by IWEs
                            for (i in 1:(length(SSE_intervals)-1)){
                                private$interval_evol(interval_length = SSE_intervals[i], dt = dt)
                                branchEvolInfo$infoIWE[[i]] <- private$singleStr[[branchEvolInfo$islands[i]]]$IWE_evol(testing = testing)
                            }
                                        # Sequence evolution along the last branch interval
                            private$interval_evol(interval_length = SSE_intervals[length(SSE_intervals)], dt = dt)
                        }
                    }
                    if(testing){
                        return(branchEvolInfo)
                    }
                }
                )
              )

#' Split Newick Tree
#'
#' This function splits a tree in Newick format into subtrees.
#' Each subtree unit can be a leaf or another subtree.
#'
#' @param tree A string in Newick format (with or without the ending semicolon).
#' @return A data frame, rows corresponding to each subtree unit and columns
#' $unit: chr (unit name), $brlen: num (unit branch length).
#' @noRd
split_newick <- function(tree) {
  if (!is.character(tree)) {stop("Input tree must be a character string.")}
  s <- strsplit(tree, '')[[1]]
  if(s[1]!="(") stop("newick tree does not begin with '('.")
  # withinU is the index of the within unit separator character ':'
  # it separates unit (leaf/subtree) name and branch length
  withinU <- c()
  # betweenU are the indexes of the between unit separator characters ',', '(' or ')'
  # it separates the units of the current subtree
  betweenU <- 1    # Initialized with unit start '(' index
  parcount <- 1    ## how many more '(' than ')' were there?##################
  for(i in 2:length(s)) { # From the first position after unit start
    if(s[i] == '(') {
      # parcount takes integer values reflecting the depth of the current unit
      # e.g. first level subtree: 1, subtree within first subtree: 2 ..
      # and 0 when the the end of the current subtree is reached ')'
      parcount <- parcount + 1
    } else if(s[i] == ')') {
      parcount <- parcount - 1
      if(parcount == 0) {
        betweenU <- c(betweenU, i)
      }
    } else if(s[i] == ':') {
      if(parcount == 1) { # Saves only withinU separator indexes ':' in first level subtree
        withinU <- c(withinU, i)
      }
    } else if(s[i] == ',') {
      if(parcount == 1) { # Saves only betweenU separator indexes affecting first level subtree
        betweenU <- c(betweenU, i)
      }
    }
  }
  # Create the data frame to save unit names and branch lengths of the current level subrtee
  L <- data.frame("unit"=rep("", length(withinU)),
                  "brlen"=rep(NA, length(withinU)),
                  stringsAsFactors = FALSE)
  for(j in 1:length(withinU)) { # For each unit (tip or subtree)
    #Uid_start is unit name start index
    Uid_start <- betweenU[j]+1 # From next character from unit start '(' look to the right
    while(s[Uid_start] == ' ') {
      Uid_start <- Uid_start+1
    }
    #Uid_start is unit name start index
    Uid_end <- withinU[j]-1 # From within unit separator ':' look to the left
    while(s[Uid_end] == ' ') {
      Uid_end <- Uid_end-1
    }
    #brlen_start is branch length start index
    brlen_start <- withinU[j]+1 # From within unit separator ':' look to the right
    while(s[brlen_start] == ' ') {
      brlen_start <- brlen_start+1
    }
    #brlen_start is branch length end index
    brlen_end <- betweenU[j+1]-1 # From next character from unit separator ', or )' look to the left
    while(s[brlen_end] == ' ') {
      brlen_end <- brlen_end-1
    }
    L[j, 1] <- paste0(s[Uid_start:Uid_end], collapse="")
    L[j, 2] <- as.numeric(paste0(s[brlen_start:brlen_end], collapse=""))
  }
  return(L)
}

#' @title treeMultiRegionSimulator
#' @importFrom R6 R6Class
#'
#' @description
#' an R6 class representing the methylation state of GpGs in different genomic
#' structures in the nodes of a tree.
#'
#' The whole CpG sequence is an object of class combiStructureGenerator.
#' Each genomic structure in it is contained in an object of class singleStructureGenerator.
#'
treeMultiRegionSimulator <- R6Class("treeMultiRegionSimulator",
                             public = list(
                               #' @field Branch Public attribute: List containing objects of class combiStructureGenerator
                               Branch=NULL,
                               #' @field branchLength Public attribute: Vector with the corresponding branch lengths of each $Branch element
                               branchLength=NULL,

                               #' @description
                               #' Simulate CpG dinucleotide methylation state evolution along a tree.
                               #' The function splits a given tree and simulates evolution along its
                               #' branches. It recursively simulates evolution in all of the subtrees in the given tree
                               #' until the tree leafs
                               #'
                               #' @param Tree String. Tree in Newick format. When called recursivelly it is given the corresponding subtree.
                               #' @param dt Length of SSE time steps.
                               #' @param testing Default FALSE. TRUE for testing purposes.
                               #' @param parent_index Default 1. When called recursivelly it is given the corresponding parent branch index.
                               #'
                               #' @return NULL
                               treeEvol = function(Tree, dt=0.01, parent_index=1, testing=FALSE) {
                                 tl <- split_newick(Tree)    ## list of subtrees and lengths of branches leading to them
                                 if(testing) {
                                   message("List of subtrees and lengths of branches leading to them")
                                   message(tl)
                                   message("Simulation on first branch")
                                 }
                                 for(i in 1:nrow(tl)) {
                                   ni <- length(self$Branch) + 1  ## "new index", that is index of the new branch
                                   self$Branch[[ni]] <<- self$Branch[[parent_index]]$copy()
                                   self$Branch[[ni]]$set_IWE_events(NULL)
                                   self$Branch[[ni]]$set_offspring_index(NULL)
                                   self$Branch[[ni]]$set_own_index(ni)
                                   self$Branch[[ni]]$set_parent_index(parent_index)
                                   self$Branch[[parent_index]]$add_offspring_index(ni)
                                   s <- tl[i, 1]
                                   self$branchLength[ni] <- tl[i, "brlen"]
                                   if(substr(s, 1, 1) == "(") {   ## real subtree
                                     self$Branch[[ni]]$branch_evol(tl[i, "brlen"], dt=dt)
                                     self$treeEvol(Tree = s, dt=dt, parent_index=ni, testing=testing)
                                   } else {   ## branch
                                     self$Branch[[ni]]$set_name(s)
                                     self$Branch[[ni]]$branch_evol(tl[i, "brlen"], dt=dt)
                                   }
                                 }
                               },

                               #' @description
                               #' Create a new treeMultiRegionSimulator object.
                               #' $Branch is a list for the tree branches, its first element represents the tree root.
                               #'
                               #' Note that one of either infoStr or rootData needs to be given. Not both, not neither.
                               #'
                               #' @param rootData combiStructureGenerator object. When given, the simulation uses its parameter values.
                               #' @param tree tree
                               #' @param infoStr  A data frame containing columns 'n' for the number of sites, and 'globalState' for the favoured global methylation state.
                               #' If initial equilibrium frequencies are given the dataframe must contain 3 additional columns: 'u_eqFreq', 'p_eqFreq' and 'm_eqFreq'
                               #' @param params Default NULL. When given: data frame containing model parameters. Note that rootData is given, its parameter values are used.
                               #' @param dt length of the dt time steps for the SSE evolutionary process
                               #' @param testing Default FALSE. TRUE for testing output.
                               #'
                               #' @return A new `treeMultiRegionSimulator` object.
                               initialize = function(infoStr = NULL, rootData = NULL, tree, params = NULL, dt = 0.01, testing = FALSE) {
                                 self$Branch <- list()
                                 if(!is.null(infoStr) && is.null(rootData)){
                                   message(paste("Simulating data at root and letting it evolve along given tree: ", tree))
                                   self$Branch[[1]] <- combiStructureGenerator$new(infoStr, params = params, testing = testing)
                                 }
                                 if(!is.null(rootData)&& is.null(infoStr)){
                                   message(paste("Simulating evolution of given data at root along given tree: ", tree))
                                   self$Branch[[1]] <- rootData$copy()
                                 }
                                 self$branchLength[1] <- NULL
                                 self$Branch[[1]]$set_own_index(1)
                                 self$treeEvol(Tree=tree, dt = dt, testing=testing)
                               }
                             )
)

