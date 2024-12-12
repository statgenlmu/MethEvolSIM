library(R6)

singleStructureGenerator$set("public","get_seq_length",function(){
  return(length(private$seq))
})

singleStructureGenerator$set("public","get_myCombiStructure",function(){
  return(length(private$my_combiStructure))
})

singleStructureGenerator$set("public","get_myCombiID",function(){
  return(private$my_combiStructure$get_id())
})


#' @description
#' The function creates a transition matrix based on the old and new equilibrium frequencies passed as input parameter
#'
#' @param old_eqFreqs Vector containing three old equilibrium frequencies adding up to 1
#' @param new_eqFreqs Vector containing three new equilibrium frequencies adding up to 1
#' @param testing Default FALSE. When TRUE warnings are thrown when the output matrix is not fulfilling its mathematical properties
#' 
#' @return A Transition Matrix calculated from the old and new equilibrium frequencies
singleStructureGenerator$set("public","get_TransMatrix",function(old_eqFreqs,new_eqFreqs, testing=FALSE){
  #set the transition matrix for the equilibrium frequencies
  m <- NULL
  
  u_new <- new_eqFreqs[1]
  p_new <- new_eqFreqs[2]
  m_new <- new_eqFreqs[3]
  
  u_old <- old_eqFreqs[1]
  p_old <- old_eqFreqs[2]
  m_old <- old_eqFreqs[3]
  
  
  # Check Case 1: 1 new frequency value bigger and 2 smaller
  
  if (u_new > u_old & p_new <= p_old & m_new <= m_old) {
    IWE_case <- "Case 1. u bigger"
    m <- matrix(c(1, 0, 0,
                  (p_old-p_new)/p_old, p_new/p_old, 0,
                  (m_old-m_new)/m_old, 0, m_new/m_old),
                nrow = 3, byrow = TRUE)
    
  }
  if (p_new > p_old & u_new <= u_old & m_new <= m_old) {
    IWE_case <- "Case 1. p bigger"
    m <- matrix(c(u_new/u_old, (u_old-u_new)/u_old, 0,
                  0, 1, 0,
                  0, (m_old-m_new)/m_old, m_new/m_old),
                nrow = 3, byrow = TRUE)
  }
  if (m_new > m_old & p_new <= p_old & u_new <= u_old) {
    IWE_case <- "Case 1. m bigger"
    m <- matrix(c(u_new/u_old, 0, (u_old-u_new)/u_old,
                  0, p_new/p_old, (p_old-p_new)/p_old,
                  0, 0, 1),
                nrow = 3, byrow = TRUE)
    
  }
  # Check Case 2: 1 new frequency value smaller
  if (u_new < u_old & p_new >= p_old & m_new >= m_old) {
    IWE_case <- "Case 2. u smaller"
    m <- matrix(c(u_new/u_old, (p_new-p_old)/u_old, (m_new-m_old)/u_old,
                  0, 1, 0,
                  0, 0, 1),
                nrow = 3, byrow = TRUE)
  }
  if (p_new < p_old & u_new >= u_old & m_new >= m_old) {
    IWE_case <- "Case 2. p smaller"
    m <- matrix(c(1, 0, 0,
                  (u_new-u_old)/p_old, p_new/p_old, (m_new-m_old)/p_old,
                  0, 0, 1),
                nrow = 3, byrow = TRUE)
  }
  if (m_new < m_old & p_new >= p_old & u_new >= u_old) {
    IWE_case <- "Case 2. m smaller"
    m <- matrix(c(1, 0, 0,
                  0, 1, 0,
                  (u_new-u_old)/m_old, (p_new-p_old)/m_old, m_new/m_old),
                nrow = 3, byrow = TRUE)
  }
  
  
  
  if(testing == TRUE){
    # Validate transition matrix
    validationStates <- listTransitionMatrix_validation(list(m), "m")
    listMatrices_validationResults(validationStates)
    # Validate Markov Chain State Transition Property
    validationStates <- transPropMC_validation(old_eqFreqs, m, new_eqFreqs, listName = "transPropMC get_TransMatrix")
    transPropMC_validationResults(validationStates)
  }
  
  
  return(m)
})

#' @name RE_within
#' @title RE_within method for recombination in ssgs
#' @method RE_within
#' @description
#' The function simulates a recombination event within a SingleStructureGenerator object
#'
#' @param Y_seq A vector containing the sequence of methylation states of the recombining SingleStructureGenerator
#' @param Y_eqFreqs Vector containing three equilibrium frequencies from the recombining SingleStructureGenerator
#' @param testing Default FALSE. 
#' 
#' @return If testing = TRUE Returns a list of the sampled equilibrium frequencies and the sequence, their id, and the new sequence after recombination.
#' Additionally warnings are thrown when rate matrices are not fulfilling their mathematical properties. 
#' \describe{
#'       \item{\code{chosenEqFreqs}}{Equilibrium frequencies chosen within the RE_within method.}
#'       \item{\code{id}}{ID of the chosen Sequence. FreqsX denotes that the frequencies of the current singleStructureGenerator were chosen, FreqsY denote they were chosen from the right sequence outside of the singlesStructureGenerator.}
#'       \item{\code{chosenSeq}}{The chosen sequence corresponding to the id.}
#'       \item{\code{newSeq}}{The newly generated sequence after the recombination event between the two SingleStructureGenerator objects}
#'      } 
#' @export

singleStructureGenerator$set("public","RE_within",function(Y_seq,Y_eqFreqs,testing = FALSE){
  new_whole_seq <- NULL
  if (testing){
    # compute previous observed methylation frequencies
    old_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
    # Save previous rates
    old_Q <- private$Q
    # Initiate changedPos with NULL
    changedPos <- NULL
  }
  #get the length of incoming sequence
  n_Y <- length(Y_seq)
  # set the length of the total sequence
  n_T <- length(private$seq)
  # set the length of the staying original sequence
  n_X <- n_T - n_Y
  # choose which are the new eqfreqs
  chosen_eqFreqs <- sample(c("FreqsX", "FreqsY"), size = 1, prob = c(n_X/n_T, n_Y/n_T))
  # get the X sequence. if recombination happens at position 1, x_seq <- c()
  X_seq <- if (n_X > 0) private$seq[1:n_X] else c()
  
  #define equilibrium frequencies for X, Y_eqFreqs is given to the function
  X_eqFreqs <- private$eqFreqs
  X_seq_testing <- NULL
  Y_seq_testing <- NULL
  X_eqFreqs_testing <- NULL
  old_eqFreqs <- NULL
  new_eqFreqs <- NULL
  seq_to_update <- NULL
  #if frequencies are identical just concatenate the sequences
  if (identical(X_eqFreqs, Y_eqFreqs)) {
    new_whole_seq <- c(X_seq, Y_seq)
  }
  
  #if the recombination happens at the first position of the ssg
  else if(length(X_seq)==0){
    private$eqFreqs <- Y_eqFreqs
    new_whole_seq <- Y_seq
    # Update Qi and Q 
    private$set_Qi()
    private$set_Q()
    if(testing){
      X_seq_testing <- X_seq
      Y_seq_testing <- Y_seq
    }
  }
  else{
    # determine which frequencies got chosen
    if (chosen_eqFreqs == "FreqsX"){
      #X_eqFreqs is the new frequencies since they get applied to Y as "new" frequencies
      #Y_eqFreqs is the old frequencies since they get replaced by X_eqFreqs
      old_eqFreqs <- Y_eqFreqs
      new_eqFreqs <- X_eqFreqs
      seq_to_update <- Y_seq
    }
    
    else if (chosen_eqFreqs =="FreqsY"){
      old_eqFreqs <- X_eqFreqs
      new_eqFreqs <- Y_eqFreqs
      seq_to_update <- X_seq
    }
    
    #compute the transition matrix for equilibrium frequencies
    Sr <- self$get_transMat(old_eqFreqs = old_eqFreqs, new_eqFreqs=new_eqFreqs,info="RE_within",testing=FALSE)
    # Change data equilibrium frequencies
    private$eqFreqs <<- new_eqFreqs 
    # Update Qi and Q 
    private$set_Qi()
    if(testing){
      X_seq_testing <- X_seq
      Y_seq_testing <- Y_seq
      # Validate updated rate matrices
      validationStates <- listRateMatrix_validation(private$Qi, "Qi matrices RE_within")
      listMatrices_validationResults(validationStates)
      
      validationStates <- listRateMatrix_validation(private$Q[[1]], "Q matrices RE_within R1")
      listMatrices_validationResults(validationStates)
      
      validationStates <- listRateMatrix_validation(private$Q[[2]], "Q matrices RE_within R2")
      listMatrices_validationResults(validationStates)
      
      validationStates <- listRateMatrix_validation(private$Q[[3]], "Q matrices RE_within R3")
      listMatrices_validationResults(validationStates)
      
      # Validate methylation equilibrium triple (extra control at RE_within)
      validationStates <- listFreqVector_validation(list(old_eqFreqs, new_eqFreqs), "RE_within control old_eqFreqs and new_eqFreqs")
      listFreqVector_validationResults(validationStates)
      
    }
    # Sample Y_seq according to transition probablities
    new_seq <- rep(0, length(seq_to_update))
    for(i in 1:length(new_seq)) {
      new_seq[i] <- sample(1:3, size=1, prob=as.vector(Sr[seq_to_update[i],]))
    }
    if (chosen_eqFreqs=="FreqsX"){
      new_whole_seq <- c(X_seq, new_seq)
    }
    else{
      new_whole_seq <- c(new_seq, Y_seq)
    }
  }
  
  if(any(private$seq != new_whole_seq)){
    changedPos <- which(private$seq != new_whole_seq)
    private$seq <- new_whole_seq
    for (i in changedPos){
      private$update_neighbSt(i)
      #delete for loop (sara)
      # Update $ratetree
      for(j in max(i-1, 1):min(i+1, length(private$seq))) {
        private$update_ratetree(j, abs(private$Q[[private$siteR[j]]][[private$neighbSt[j]]][private$seq[j],private$seq[j]]))
      }
      #update instead of for loop beaause of changes to MethEvolSIM  
      #  if(i==1){
      #    private$update_ratetree_betweenStr(prevStr=TRUE)
      #  }
      #  if(i==length(private$seq)){
      #    private$update_ratetree_betweenStr(nextStr=TRUE)
      #    
      #  }
    }
  }
  #self$initialize_ratetree()
  
  #Compute the new_obsFreqs after the RE event
  new_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
  if(testing){
    if(chosen_eqFreqs=="FreqsX"){
      return(list(chosenEqFreqs = X_eqFreqs, id="FreqsX", chosenSeq = X_seq_testing, newSeq = new_whole_seq))
    }
    else{
      return(list(chosenEqFreqs = Y_eqFreqs,id= "FreqsY",chosenSeq = Y_seq_testing, newSeq = new_whole_seq))
    }
  }
  
  
})

#' @description
#' The function simulates a recombination event within a CombiStructureGenerator object
#' @name RE_complete
#' @title Method for executing recombination in a csg
#' @param combi_2 A CombiStructureGenerator object that recombines with the current CombiStructureGenerator
#' @param position Integer of the recombination position inside the CombiStructureGenerator Object
#' @param testing Default FALSE. 
#' 
#' @return If testing = TRUE Returns a list containing the following items:
#' \describe{
#'       \item{\code{combi_2}}{The recombining CombiStructureGenerator object.}
#'       \item{\code{combi_1}}{The current CombiStructureGenerator object.}
#'       \item{\code{recomb_point}}{The recombination position also referred as recombination point.}
#'       \item{\code{ssg_index}}{The index of the SingleStructureGenerator object where the recombination is happening.}
#'       \item{\code{ssg_position_within}}{The relative recombination position within the SingleStructureGenerator object where the recombination happens.}
#'       \item{\code{recombEqFreqs}}{Equilibrium frequencies chosen within the SingleStructureGenerator object where the recombination happens.}
#'       \item{\code{new_ssgs}}{A list of all SingleStructureGenerator objects contained in the current CombiSturctureGenerator object after the recombination event.}
#'      } 


combiStructureGenerator$set("public","RE_complete",function(combi_2, position,testing=FALSE){
  #initializing variables to keep track of the length and the target single structure
  cumulative_length <- 0
  target_singleStr_index <- NULL
  target_position_in_singleStr <- NULL
  testing_info <- list()
  if(testing){
    testing_info <- append(testing_info,list(combi_2 = combi_2, combi_1 = self, recomb_point = position))
  }
  
  if (self$get_singleStr_number()==1){
    Y_seq <- combi_2$get_singleStr(1)$get_seq()[position:length(combi_2$get_singleStr(1)$get_seq())]
    Y_eqFreqs <- combi_2$get_singleStr(1)$get_eqFreqs()
    private$singleStr[[1]]$RE_within(Y_seq = Y_seq, Y_eqFreqs = Y_eqFreqs)
    
  }
  
  else{
    #loop through each stinglestructure in the combistructure
    for(i in 1:self$get_singleStr_number()){
      singleStr_length <- private$singleStr[[i]]$get_seq_length()
      cumulative_length <- cumulative_length + singleStr_length
      
      #check if the recombination position falls within the current singlestr
      if(cumulative_length >= position){
        target_singleStr_index <- i
        target_position_in_singleStr <- position-(cumulative_length-singleStr_length)
        break 
      }
    }
    combi_2_target_ssg <- combi_2$get_singleStr(target_singleStr_index)
    Y_seq <- combi_2_target_ssg$get_seq()[(target_position_in_singleStr):length(combi_2_target_ssg$get_seq())]
    Y_eqFreqs <- combi_2_target_ssg$get_eqFreqs()
    out <- private$singleStr[[target_singleStr_index]]$RE_within(Y_seq = Y_seq, Y_eqFreqs = Y_eqFreqs,testing=testing)
    
    if(testing){
      testing_info <- append(testing_info, list(ssg_index = target_singleStr_index, ssg_position_within = target_position_in_singleStr,recombEqFreqs = out$chosenEqFreqs))
    }
    
    #when the recombination happens in the last ssg then this step is not needed since RE_within did everything
    if(target_singleStr_index < combi_2$get_singleStr_number()){
      for(j in (target_singleStr_index+1):combi_2$get_singleStr_number()){
        private$singleStr[[j]] <- combi_2$copy()$get_singleStr(j)
        private$singleStr[[j]]$set_myCombiStructure(self)
      }
    }
    
  }
  
  if(testing){
    testing_info <- append(testing_info,list(new_ssgs=private$singleStr))
    return(testing_info)
  }
  
})



recombMultiRegionSimulator <- R6Class("recombMultiRegionSimulator", public = list(
  nodes = NULL,
  edges = NULL,
  node_list = NULL,
  tip_list = NULL,
  testing = FALSE,
  infoStr = NULL,
  params = NULL,
  testing_list = list(order = list(), combis = list()),
  
  initialize = function(nodes, edges, infoStr = NULL, params = NULL, dt = 0.01, testing){
    #Ensure nodes and edges are data frames and contain the necessary columns
    stopifnot(is.data.frame(nodes), is.data.frame(edges))
    stopifnot("time" %in% colnames(nodes), "id" %in% colnames(nodes))
    stopifnot("parent" %in% colnames(edges), "child" %in% colnames(edges), "left" %in% colnames(edges), "right" %in% colnames(edges))
    
    self$params <- params
    self$infoStr <- infoStr
    self$testing <- testing
    self$nodes <- nodes
    self$edges <- edges
    self$nodes$id <- self$nodes$id + 1
    self$edges$parent <- self$edges$parent + 1
    self$edges$child <- self$edges$child + 1
    self$classify_and_order_events()
    

    self$node_list <- rep(list(NULL),nrow(nodes))
    self$evolve()
  },
  
  classify_and_order_events = function(){
    self$nodes <- self$nodes[order(self$nodes$time,decreasing = TRUE),]
  },
  
  evolve = function(dt=0.01,testing=self$testing){
    testing_combis <- list(combi_1 = list())
    #a list containing all the existing nodes that couldnt be recombined
    pending_recombine <- rep(list(NULL),nrow(self$nodes))
    test <- 0
    for(i in 1:nrow(self$nodes)){
      #if node recombination is pending skip it
      if(is.null(pending_recombine[[self$nodes$id[[i]]]])& self$nodes$time[[i]]>0){
        parent_id <- self$nodes$id[[i]]
        parent_time <- self$nodes$time[[i]]
        if(nrow(self$edges[self$edges$child == parent_id,])==0){
          #new root identified
          test <- test+1
          #message("New root ID: ", parent_id-1, " Time: ", self$nodes[self$nodes$id==parent_id,"time"], "Root number: ", test)
          self$node_list[[parent_id]] <- list(combi = combiStructureGenerator$new(infoStr = self$infoStr, params = self$params, testing = testing), right = NULL)
          if(testing){
            self$testing_list$order[[parent_id]] <- 0
          }
          
        }
        #for the recombining node that doesnt get further evolved we dont need to follow the tree further
        outgoing_edges <- self$edges[self$edges$parent == parent_id,]
        for(edge_row_index in 1:nrow(outgoing_edges)){
          child_id <- outgoing_edges$child[[edge_row_index]]
          child_time <- self$nodes[self$nodes$id == outgoing_edges$child[[edge_row_index]],"time"]
          #check wether the child was already processed, this happens when two nodes point to the same child node. 
          #this would be a problem when that child node also is recombining because we then execute recombination twice
          if(is.null(self$node_list[[child_id]])){
            self$node_list[[child_id]] <- list(combi = self$node_list[[parent_id]]$combi$copy(), right = outgoing_edges$right[[edge_row_index]])
            self$node_list[[child_id]]$combi$branch_evol(branch_length = parent_time - child_time, dt = 0.001)

            #if recombination happens
            if(self$nodes[self$nodes$id == child_id,"flags"]==131072 & is.null(pending_recombine[[parent_id]])){
              #check which the second recombining node is in the table (lower or upper <- +/-1)
              rec_node_id <- NULL
              #alternative to looking in the dataframe is looking in pending_recombine
              if(self$nodes[self$nodes$id==child_id+1,"time"] == child_time){
                rec_node_id <- child_id+1
              }
              else if(self$nodes[self$nodes$id==child_id-1,"time"] == child_time){
                rec_node_id <- child_id-1
              }
              else{
                warning("Warning: No recombination partner found")
              }
              
              if(!is.null(self$node_list[[rec_node_id]])){
                #recombination happens
                #i always choose the right of the first combi as recombination point
                rec_position <- min(outgoing_edges$right[[edge_row_index]],self$node_list[[rec_node_id]]$right)
                out <- self$node_list[[child_id]]$combi$RE_complete(combi_2=self$node_list[[rec_node_id]]$combi,position=rec_position, testing = testing)
                if(testing){
                  self$testing_list$order[[child_id]] <- 131072
                  self$testing_list$order[[rec_node_id]] <- 131072
                  testing_combis$combi_1[[length(testing_combis$combi_1)+1]] <- out$combi_1

                }
                
                #message("recombination between: ",rec_node_id-1, " and ", child_id-1)
                pending_recombine[rec_node_id] <- list(NULL)
                #when recombination happens delete from pending_recombine and add to processed nodes
              }
              else{
                #recombination cannot happen rn since recombining node is not evolved yet
                pending_recombine[[child_id]] <- child_id
              }
              
            }
            else{
              if(testing & !(self$nodes[self$nodes$id == child_id,"flags"]==131072) || !(is.null(pending_recombine[[parent_id]]))){
                self$testing_list$order[[child_id]] <- 0
                
              }
            }
          }
          
          
        }
      } 
      else{
        if(self$nodes$time[[i]]==0){
          self$tip_list[[self$nodes$id[[i]]]] <- self$node_list[[self$nodes$id[[i]]]]
        }  
        #message("Skip Case: pending_recombine[[i]] = ",pending_recombine[[i]]," self$nodes$time[[i]] = ", self$nodes$time[[i]], " i = ",i," self$nodes$id[[i]] = ",self$nodes$id[[i]])
        if(testing){
          self$testing_list$order[[self$nodes$id[[i]]]] <- 1
          self$testing_list$combis <- testing_combis
        }
      }
    }
  }
  
))


#example to run the simulation. a list of the simulated ARG tips sequences gets printed. 
'
set.seed(12212333)
nodes <- read.csv("path/to/nodes")
edges <- read.csv("path/to/edges")
prams <- data.frame(  alpha_pI= 0.1,
                      beta_pI= 1,
                      alpha_mI= 0.1,
                      beta_mI= 0.5,
                      alpha_pNI= 0.1,
                      beta_pNI= 1,
                      alpha_mNI= 0.5,
                      beta_mNI= 0.1,
                      alpha_Ri= 0.1,
                      iota= 0.3,
                      mu= 0.26
                      )

infoStr <- data.frame(
  n = c(200,200,200,200,200),
  globalState = c("M","U","M","U","M")  
)

recombination_handler <- recombMultiRegionSimulator$new(
  nodes = nodes,
  edges = edges,
  infoStr = infoStr,
  params = prams,
  dt = 0.01,
  testing = TRUE
)
tip_list <- recombination_handler$tip_list
for(i in 1:number_of_tis){
  print(tip_list[[i]]$combi$get_singleStr(4)$get_seq())
}
print("-------")
for(i in 1:number_of_tips){
  print(tip_list[[i]]$combi$get_singleStr(2)$get_seq())
}
'



