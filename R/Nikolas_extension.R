library(R6)
library(devtools)


singleStructureGenerator$set("public","get_eqMatrix",function(old_eqFreqs,new_eqFreqs, testing=FALSE){
  #set the transition matrix for the equilibrium frequencies
  Mk <- NULL
  
  u_new <- new_eqFreqs[1]
  p_new <- new_eqFreqs[2]
  m_new <- new_eqFreqs[3]
  
  u_old <- old_eqFreqs[1]
  p_old <- old_eqFreqs[2]
  m_old <- old_eqFreqs[3]

  
  # Check Case 1: 1 new frequency value bigger and 2 smaller
  
  if (u_new > u_old & p_new <= p_old & m_new <= m_old) {
    IWE_case <- "Case 1. u bigger"
    Mk <- matrix(c(1, 0, 0,
                   (p_old-p_new)/p_old, p_new/p_old, 0,
                   (m_old-m_new)/m_old, 0, m_new/m_old),
                 nrow = 3, byrow = TRUE)
    
  }
  if (p_new > p_old & u_new <= u_old & m_new <= m_old) {
    IWE_case <- "Case 1. p bigger"
    Mk <- matrix(c(u_new/u_old, (u_old-u_new)/u_old, 0,
                   0, 1, 0,
                   0, (m_old-m_new)/m_old, m_new/m_old),
                 nrow = 3, byrow = TRUE)
  }
  if (m_new > m_old & p_new <= p_old & u_new <= u_old) {
    IWE_case <- "Case 1. m bigger"
    Mk <- matrix(c(u_new/u_old, 0, (u_old-u_new)/u_old,
                   0, p_new/p_old, (p_old-p_new)/p_old,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
    
  }
  # Check Case 2: 1 new frequency value smaller
  if (u_new < u_old & p_new >= p_old & m_new >= m_old) {
    IWE_case <- "Case 2. u smaller"
    Mk <- matrix(c(u_new/u_old, (p_new-p_old)/u_old, (m_new-m_old)/u_old,
                   0, 1, 0,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
  }
  if (p_new < p_old & u_new >= u_old & m_new >= m_old) {
    IWE_case <- "Case 2. p smaller"
    Mk <- matrix(c(1, 0, 0,
                   (u_new-u_old)/p_old, p_new/p_old, (m_new-m_old)/p_old,
                   0, 0, 1),
                 nrow = 3, byrow = TRUE)
  }
  if (m_new < m_old & p_new >= p_old & u_new >= u_old) {
    IWE_case <- "Case 2. m smaller"
    Mk <- matrix(c(1, 0, 0,
                   0, 1, 0,
                   (u_new-u_old)/m_old, (p_new-p_old)/m_old, m_new/m_old),
                 nrow = 3, byrow = TRUE)
  }
  return(Mk)
  
  
})



singleStructureGenerator$set("public","RE_within",function(Y_seq,Y_eqFreqs,testing = FALSE){
  print(private$seq)
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
  
  # determine which frequencies got chosen
  if (chosen_eqFreqs == "FreqsX"){
    #X_eqFreqs is the new frequencies since they get applied to Y as "new" frequencies
    #Y_eqFreqs is the old frequencies since they get replaced by X_eqFreqs
    
    if (identical(X_eqFreqs, Y_eqFreqs)) {
      eqFreqsChange = F
      new_whole_seq <- c(X_seq, Y_seq)
    } else{
      eqFreqsChange = T
      #compute the transition matrix for equilibrium frequencies
      Mk <- self$get_eqMatrix(old_eqFreqs = Y_eqFreqs, new_eqFreqs=X_eqFreqs)
      
      # Change data equilibrium frequencies
      private$eqFreqs <<- X_eqFreqs
      # Update Qi and Q
      private$set_Qi()
      private$set_Q()
      
      # Sample Y_seq according to transition probablities
      new_Y_seq <- rep(0, length(Y_seq))
      for(i in 1:length(new_Y_seq)) {
        new_Y_seq[i] <- sample(1:3, size=1, prob=as.vector(Mk[Y_seq[i],]))
      }
      new_whole_seq <- c(X_seq, new_Y_seq)
    }
    
    
    # Update $seq and $neighbSt
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
  }
  else if (chosen_eqFreqs =="FreqsY"){
    #the new equilibrium frequencies here are the Y frequencies
    #the old equilibrium frequencies here are the X frequencies
    
    if (identical(X_eqFreqs, Y_eqFreqs)) {
      eqFreqsChange = F
      new_whole_seq <- c(X_seq, Y_seq)
    } 
    else{
      eqFreqsChange = T
      #compute the transition matrix for equilibrium frequencies
      Mk <- self$get_eqMatrix(old_eqFreqs = X_eqFreqs, new_eqFreqs=Y_eqFreqs)
      
      # Change data equilibrium frequencies
      private$eqFreqs <<- Y_eqFreqs
      # Update Qi and Q
      private$set_Qi()
      private$set_Q()
      
      # Sample X_seq according to transition probablities
      new_X_seq <- rep(0, length(X_seq))
      for(i in 1:length(new_X_seq)) {
        new_X_seq[i] <- sample(1:3, size=1, prob=as.vector(Mk[X_seq[i],]))
      }
      new_whole_seq <- c(new_X_seq, Y_seq)
    }
    
    
    # Update $seq and $neighbSt
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
  }
  
  #Compute the new_obsFreqs after the RE event
  new_obsFreqs <- c(sum(private$seq==1), sum(private$seq==2), sum(private$seq==3))/length(private$seq)
  print(new_obsFreqs)
  print(private$seq)
})


#for testing
#ssg <- singleStructureGenerator$new(globalState = "U", n = 40,eqFreqs = c(1,0,0))
#print(ssg$RE_within(Y_seq=c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),Y_eqFreqs=c(0.7,0.1,0.2)))



