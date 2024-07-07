library(R6)
library(devtools)

singleStructureGenerator$set("public","get_length",function(){
  return(length(private$seq))
})

singleStructureGenerator$set("public","get_TransMatrix",function(old_eqFreqs,new_eqFreqs, testing=FALSE){
  #set the transition matrix for the equilibrium frequencies
  Sr <- NULL
  
  u_new <- new_eqFreqs[1]
  p_new <- new_eqFreqs[2]
  m_new <- new_eqFreqs[3]
  
  u_old <- old_eqFreqs[1]
  p_old <- old_eqFreqs[2]
  m_old <- old_eqFreqs[3]
  
  
  # Check Case 1: 1 new frequency value bigger and 2 smaller
  
  if (u_new > u_old & p_new <= p_old & m_new <= m_old) {
    IWE_case <- "Case 1. u bigger"
    Sr <- matrix(c(1, 0, 0,
                   (p_old-p_new)/p_old, p_new/p_old, 0,
                   (m_old-m_new)/m_old, 0, m_new/m_old),
                 nrow = 3, byrow = TRUE)
    
  }
  if (p_new > p_old & u_new <= u_old & m_new <= m_old) {
    IWE_case <- "Case 1. p bigger"
    Sr <- matrix(c(u_new/u_old, (u_old-u_new)/u_old, 0,
                   0, 1, 0,
                   0, (m_old-m_new)/m_old, m_new/m_old),
                 nrow = 3, byrow = TRUE)
  }
  if (m_new > m_old & p_new <= p_old & u_new <= u_old) {
    IWE_case <- "Case 1. m bigger"
    Sr <- matrix(c(u_new/u_old, 0, (u_old-u_new)/u_old,
                   0, p_new/p_old, (p_old-p_new)/p_old,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
    
  }
  # Check Case 2: 1 new frequency value smaller
  if (u_new < u_old & p_new >= p_old & m_new >= m_old) {
    IWE_case <- "Case 2. u smaller"
    Sr <- matrix(c(u_new/u_old, (p_new-p_old)/u_old, (m_new-m_old)/u_old,
                   0, 1, 0,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
  }
  if (p_new < p_old & u_new >= u_old & m_new >= m_old) {
    IWE_case <- "Case 2. p smaller"
    Sr <- matrix(c(1, 0, 0,
                   (u_new-u_old)/p_old, p_new/p_old, (m_new-m_old)/p_old,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
  }
  if (m_new < m_old & p_new >= p_old & u_new >= u_old) {
    IWE_case <- "Case 2. m smaller"
    Sr <- matrix(c(1, 0, 0,
                   0, 1, 0,
                   (u_new-u_old)/m_old, (p_new-p_old)/m_old, m_new/m_old),
                 nrow = 3, byrow = TRUE)
  }
  
  
  
  if(testing == TRUE){
    # Validate transition matrix
    validationStates <- listTransitionMatrix_validation(list(Sr), "Sr")
    listMatrices_validationResults(validationStates)
    # Validate Markov Chain State Transition Property
    validationStates <- transPropMC_validation(old_eqFreqs, Sr, new_eqFreqs, listName = "transPropMC get_SrMatrix")
    transPropMC_validationResults(validationStates)
  }
  
  
  return(Sr)
})



singleStructureGenerator$set("public","RE_within",function(Y_seq,Y_eqFreqs,testing = FALSE){
  #NEED TO COVER THE CASE WHERE THE YSEQ IS EQUAL OR GREATER THAN THE X SEQ
  
  new_whole_seq <- NULL
  if (testing){
    # compute previous observed methylation frequencies
    old_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
    # Save previous rates
    old_Q <- private$Q
    # Initiate changedPos with NULL
    changedPos <- NULL
    
    
    
    
    #print("Original sequence:")
    #print(private$seq)
  }
  #get the length of incoming sequence
  n_Y <- length(Y_seq)
  # set the length of the total sequence
  n_T <- length(private$seq)
  # set the length of the staying original sequence
  n_X <- n_T - n_Y
  # choose which are the new eqfreqs
  chosen_eqFreqs <- sample(c("FreqsX", "FreqsY"), size = 1, prob = c(n_X/n_T, n_Y/n_T))
  # get the X sequence
  X_seq <- private$seq[1:n_X]
  
  #define equilibrium frequencies for X, Y_eqFreqs is given to the function
  X_eqFreqs <- private$eqFreqs
  
  old_eqFreqs <- NULL
  new_eqFreqs <- NULL
  seq_to_update <- NULL
  
  #if frequencies are identical just concatenate the sequences
  if (identical(X_eqFreqs, Y_eqFreqs)) {
    eqFreqsChange = F
    new_whole_seq <- c(X_seq, Y_seq)
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
      
    eqFreqsChange = T
    #compute the transition matrix for equilibrium frequencies
    Sr <- self$get_TransMatrix(old_eqFreqs = old_eqFreqs, new_eqFreqs=new_eqFreqs)
    # Change data equilibrium frequencies
    private$eqFreqs <<- new_eqFreqs
    # Update Qi and Q
    private$set_Qi()
    private$set_Q()
      
    if(testing){
      # Validate updated rate matrices
      validationStates <- listRateMatrix_validation(private$Qi, "Qi matrices RE_within")
      listMatrices_validationResults(validationStates)
    
      validationStates <- listRateMatrix_validation(private$Q[[1]], "Q matrices RE_within R1")
      listMatrices_validationResults(validationStates)
        
      validationStates <- listRateMatrix_validation(private$Q[[2]], "Q matrices RE_within R2")
      listMatrices_validationResults(validationStates)
        
      validationStates <- listRateMatrix_validation(private$Q[[3]], "Q matrices RE_within R3")
      listMatrices_validationResults(validationStates)

      # Validate methylation equilibrium tripple (extra control at RE_within)
      validationStates <- listFreqVector_validation(list(old_eqFreqs, new_eqFreqs), "RE_within control old_eqFreqs and new_eqFreqs")
      listFreqVector_validationResults(validationStates)
      
    }
    print("1")
    print(Sr)
    print(chosen_eqFreqs)
    print(X_seq)
    print(Y_seq)
    # Sample Y_seq according to transition probablities
    new_seq <- rep(0, length(seq_to_update))
    for(i in 1:length(new_seq)) {
      print(seq_to_update)
      print(seq_to_update[i])
      print(as.vector(Sr[seq_to_update[i],]))
      new_seq[i] <- sample(1:3, size=1, prob=as.vector(Sr[seq_to_update[i],]))
    }
    if (chosen_eqFreqs=="FreqsX"){
      new_whole_seq <- c(X_seq, new_seq)
    }
    else{
      new_whole_seq <- c(new_seq, Y_seq)
    }
  }
  print("2")
  if(any(private$seq != new_whole_seq)){
    changedPos <- which(private$seq != new_whole_seq)
    private$seq <- new_whole_seq
    for (i in changedPos){
      private$update_neighbSt(i)
      # Update $ratetree
      for(j in max(i-1, 1):min(i+1, length(private$seq))) {
        private$update_ratetree(j, abs(private$Q[[private$siteR[j]]][[private$neighbSt[j]]][private$seq[j],private$seq[j]]))
      }
    }
  }
  print("3")

  #Compute the new_obsFreqs after the RE event
  new_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
  if(testing){
    if(chosen_eqFreqs=="FreqsX"){
      return(list(chosenEqFreqs = X_eqFreqs, id="FreqsX"))
    }
    else{
      return(list(chosenEqFreqs = Y_eqFreqs,id= "FreqsY"))
    }
  }
  print("3")


})

combiStructureGenerator$set("public","RE_complete",function(combi_2, position){
  #initializing variables to keep track of the length and the target single structure
  cumulative_length <- 0
  target_singleStr_index <- NULL
  target_position_in_singleStr <- NULL
  
  cumulative_length_2 <- 0
  target_singleStr_index_2 <- NULL
  target_position_in_singleStr_2 <- NULL
  
  #loop through each stinglestructure in the combistructure
  for(i in 1:self$get_singleStr_number()){
    singleStr_length <- private$singleStr[[i]]$get_length()
    cumulative_length <- cumulative_length + singleStr_length
    
    #check if the recombination position falls within the current singlestr
    if(cumulative_length >= position){
      target_singleStr_index <- i
      target_position_in_singleStr <- position-(cumulative_length-singleStr_length)
      break
    }
  }
  
  for(i in 1:combi_2$get_singleStr_number()){
    singleStr_length_2 <- combi_2$get_singleStr(i)$get_length()
    cumulative_length_2 <- cumulative_length_2 + singleStr_length_2
    
    if(cumulative_length_2 >= position){
      target_singleStr_index_2 <- i
      target_position_in_singleStr_2 <- position-(cumulative_length_2-singleStr_length_2)
      break
    }
  }
  
  combi_2_target_ssg <- combi_2$get_singleStr(target_singleStr_index_2)
  recomb_Yseq <- combi_2_target_ssg$get_seq()[target_position_in_singleStr_2:length(combi_2_target_ssg$get_seq())]
  recomb_YeqFreqs <- combi_2_target_ssg$get_eqFreqs()
  print("Test:")
  print(recomb_Yseq)
  private$singleStr[[target_singleStr_index]]$RE_within(Y_seq = recomb_Yseq, Y_eqFreqs = recomb_YeqFreqs)
  
  for(j in (target_singleStr_index+1):combi_2$get_singleStr_number()){
    private$singleStr[[j]]<<- combi_2$get_singleStr(j)$clone()
  }
  
  return(private$singleStr)

})









#for testing
#ssg <- singleStructureGenerator$new(globalState = "U", n = 40,eqFreqs = c(1,0,0))
#print(ssg$RE_within(Y_seq=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),Y_eqFreqs=c(0.7,0.1,0.2)))
infoStr <- data.frame(n = c(12, 13, 13),
                      globalState = c("M", "U", "M"))
csg <- combiStructureGenerator$new(infoStr = infoStr)

infoStr_2 <- data.frame(n= c(5,5,3),globalState=c("M","M","U"))
csg_2 <- combiStructureGenerator$new(infoStr = infoStr_2)
print(csg_2$get_singleStr_number())
for(i in 1:csg_2$get_singleStr_number()){
  print(csg_2$get_singleStr(i)$get_seq())
}


csg$RE_complete(combi_2 = csg_2,position = 8)








