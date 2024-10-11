## next steps and todos:
##
## Both of us should check whether this is really the CFTP algorithm and understand why the CFTP algorithm
## leads to a sequence sampled from the equilibrium distribution of the Castillo-Vicente--Metzler model.
##
## Testing: Unit tests should check whether applications of events as above to two sequences a and b
## of which a <= b is fulfilled ("a <= b" means that "a < b" or "a = b"), always lead to sequences a' and b'
## with a' <= b'.
##
## More time-consuming function tests could check with certain summary statistics, e.g. combinations of
## neighboring pairs or triples within CpGs, in non-CpG or at boundaries have the same frequencies in initial
## cftp-generated states as in states that evolved for a while after stating in such a state.
## In particular: is state p (partial methylation) underrepresented?
##
## Integrate contents of this file into class defintions in file multiRegion_SIM.R;
## check for functions that are currently public whether they can be private
##
## Extend CFTP by allowing also IWEs.
##
## Perhaps later: look for possibilities to improve structure of implementation (refactoring) and efficiency


################################################################################
################################################################################


## @description
## Public Method. 
## 
## @param siteR default NULL. Numerical value encoding for the sites rate of independent SSE (1, 2 or 3)
## @param oldSt default NULL. Numerical value encoding for the sites old methylation state (1, 2 or 3)
## @param newSt default NULL. Numerical value encoding for the sites new methylation state (1, 2 or 3)
## 
## @return With NULL arguments, the list of SSEi rate matrices. With non NULL arguments, the corresponding rate of change.
singleStructureGenerator$set("public", "get_Qi", function(siteR = NULL, oldSt = NULL, newSt = NULL){
  if(is.null(siteR) && is.null(oldSt) && is.null(newSt)){
    private$Qi
  } else{
    private$Qi[[siteR]][oldSt,newSt]
  }
})

##singleStructureGenerator$set("public", "get_iota", function(){
##  private$iota
##})

## @description
## Public Method.
## 
## @param index. Integer index value for the CpG position within the singleStr instance
## 
## @return decoded methylation state ($seq) of left neighbor (1, 2 or 3 for unmethylated, partially methylated or methylated)
singleStructureGenerator$set("public", "get_seqSt_leftneighb", function(index){
  if (is.null(index)) stop("Argument 'index' for CpG position is missing, with no default")
  if(!(is.numeric(index) && length(index) == 1 && index == floor(index))){
    stop("Argument 'index' must be one integer index value")
  }
  if(!(index >= 1 && index <= length(private$seq))) {
    stop("Argument 'index' not within sequence length")
  }
  # Extract row number (corresponding to $seq state of left neighbor) from the positions neighbSt
  which(private$mapNeighbSt_matrix == self$get_neighbSt(index), arr.ind = TRUE)[1]
  #floor((self$get_neighbSt(i)-1) / 3) +1
})

## @description
## Public Method.
## 
## @param index. Integer index value for the CpG position within the singleStr instance
## 
## @return decoded methylation state ($seq) of right neighbor (1, 2 or 3 for unmethylated, partially methylated or methylated)
singleStructureGenerator$set("public", "get_seqSt_rightneighb", function(index){
  if (is.null(index)) stop("Argument 'index' for CpG position is missing, with no default")
  if(!(is.numeric(index) && length(index) == 1 && index == floor(index))){
    stop("Argument 'index' must be one integer index value")
  }
  if(!(index >= 1 && index <= length(private$seq))) {
    stop("Argument 'index' not within sequence length")
  }
  which(private$mapNeighbSt_matrix == self$get_neighbSt(index), arr.ind = TRUE)[2]
  #(self$get_neighbSt(i)-1) %% 3 +1
})


## @description
## Public Method. Make a singleStructure with the same segment lengths and parameters
## as the focal one but where all states are m or u
## 
## @param state. Character value "U" or "M"
## 
## @return right neighbSt
singleStructureGenerator$set("public", "cftp_all_equal", function(state, testing = FALSE) {
  n <- length(private$seq)
  if (state == "U"){
    ### TODO: If cftp at some point accounts for IWEs, then this assigning "U" or "M" to globalState should not be done.
    private$globalState <- "U"
    private$seq <- rep(1L, n)
  } else if (state == "M"){
    private$globalState <- "M"
    private$seq <- rep(3L, n)
  } else {
    stop("Invalid 'state' value. Must be 'U' for unmethylated or 'M' for methylated")
  }
  if (testing){
    private$seq
  }
})




## @description
## Public Method. Set the methylation state of a sequence position and update the neighbor's neighbSt. It does NOT update RATETREE 
## 
## @param index. Numerical value for the index of the CpG position within the singleStr instance
## @param newSt. Numerical value encoding for the sites new methylation state (1, 2 or 3)
## @param testing default FALSE. TRUE for testing output
## 
## @return NULL when testing FALSE. Testing output when testing TRUE.
singleStructureGenerator$set("public", "set_seqSt_update_neighbSt", function(index, newSt, testing = FALSE){
  if(!(is.numeric(index) && length(index) == 1 && index == floor(index))){
    stop("'index' must be one number without decimals")
  }
  if(!(index >= 1 && index <= length(private$seq))) {
    stop("'index' not within sequence length")
  }
  if(!(newSt %in% c(1, 2, 3))){
    stop("'newSt' must be 1, 2 or 3 (for unmethylated, partially methylated or methylated)")
  }
  private$seq[index] <- newSt
  private$update_neighbSt(index)
  if (testing){
    list(seq = private$seq,
         neighbSt = private$neighbSt)
  }
})

## @description
## Public Method. Generates the events to apply for CFTP.
## 
## @param steps. Integer value >=1
## @param testing default FALSE. TRUE for testing output
## 
## @return NULL when testing FALSE. Testing output when testing TRUE.
##
## @details
## The function add steps to the existing ones. 
## If called several times the given steps need to be higher than the sum of steps generated before.
combiStructureGenerator$set("public", "cftp_event_generator", function(steps, testing = FALSE) {
  
  if(!(is.numeric(steps) && length(steps) == 1 && steps == floor(steps) && steps >= 1)){
    stop("'index' must be one number >=1 without decimals")
  }
  
  # Set CFTP_highest_rate to be the highest rate across all singleStr withing combiStr instance
  private$CFTP_highest_rate <- 0 # ensure minimum value of 0
  singleStr_n <- c()
  for(str in 1:length(private$singleStr)){
    private$CFTP_highest_rate <- max(private$CFTP_highest_rate,            
                                     max(c(unlist(private$singleStr[[str]]$get_Qi()), 
                                           (1-private$singleStr[[str]]$get_iota())/2)))
    singleStr_n[str] <- length(private$singleStr[[str]]$get_seq())
  }
  
  # Generate for each CFTP step the event and location to apply it and a threshold for acceptance/rejection
  old_steps <- length(private$CFTP_event)
  if (steps <= old_steps){stop("The given number of steps has already been generated")}
  new_steps <- steps - old_steps
  chosen_singleStr <- integer(length=new_steps)
  chosen_site <- integer(length=new_steps)
  event <- integer(length=new_steps)
  random_threshold <- numeric(length=new_steps)
  for(n in new_steps:1) {
    # For generation -n
    # Propose 1 site and what may happen to it
    chosen_singleStr[n] <- sample(1:length(private$singleStr), 1, prob=singleStr_n)
    chosen_site[n] <- sample(1:length(private$singleStr[[str]]$get_seq()), 1)
    event[n] <- sample(1:5, 1)  ## 1,2,3: go to u, p, m by SSEi ## 4,5: copy left, copy right.                        
    # Sample a threshold to accept or reject event
    random_threshold[n] <- runif(1) # numerical value between 0 and 1
  }
  private$CFTP_chosen_singleStr <- c(private$CFTP_chosen_singleStr, chosen_singleStr)
  private$CFTP_chosen_site <- c(private$CFTP_chosen_site, chosen_site)
  private$CFTP_event <- c(private$CFTP_event, event)
  private$CFTP_random <- c(private$CFTP_random, random_threshold) 
  
  if(testing){
    list(CFTP_chosen_singleStr = private$CFTP_chosen_singleStr,
         CFTP_chosen_site = private$CFTP_chosen_site,
         CFTP_event = private$CFTP_event,
         CFTP_random = private$CFTP_random)
  }
})

## @description
## Public Method. Applies the CFTP events.
## 
## @param testing default FALSE. TRUE for testing output
## 
## @return NULL when testing FALSE. Testing output when testing TRUE.
combiStructureGenerator$set("public", "cftp_apply_events", function(testing = FALSE) {
  if(length(private$CFTP_event) < 1) stop("No CFTP events generated yet.")
  if (testing){
    # Set vector to save acceptance/rejection as T or F
    event_acceptance <- logical(length = length(private$CFTP_event))
    # Set vector to save the rate of the chosen site (j) at each CFTP step (k)
    r_jk <- rep(NA, length(private$CFTP_event))
    # Set vector to save the applied case 
    applied_event <- rep(NA, length(private$CFTP_event))
  } 
  ## TODO: Update testing vectors
  for(k in length(private$CFTP_event):1) {
    ### applies the CFTP_events from -n to 1 to the combiStructure
    i <- private$CFTP_chosen_singleStr[k]
    j <- private$CFTP_chosen_site[k]
    oldSt <- private$singleStr[[i]]$get_seq()[j]
    if(private$CFTP_event[k] < 4) { # If the event is of type SSEi, set the new St as the sampled event 
      newSt <- private$CFTP_event[k]      
      if(oldSt != newSt) {
        siteR <- private$singleStr[[i]]$get_siteR(j)
        neighbSt <- private$singleStr[[i]]$get_neighbSt(j)
        r <- private$singleStr[[i]]$get_Qi(siteR = siteR, oldSt = oldSt, newSt = newSt)
        if( r/private$CFTP_highest_rate > private$CFTP_random[k] ) {
          private$singleStr[[i]]$set_seqSt_update_neighbSt(j, newSt)
          if (testing){
            event_acceptance[k] <- TRUE
            r_jk[k] <- r
            applied_event[k] <- paste("SSEi", newSt, sep = "_")
          } 
        } else if (testing){
          r_jk[k] <- r
        }
      } 
    } else {
      r <- (1-private$singleStr[[i]]$get_iota())/2 # Rc/2 for each event copy left or copy right
      if( r/private$CFTP_highest_rate > private$CFTP_random[k] ) {
         if(private$CFTP_event[k] == 4) {
           ##copy from left neighbor
           private$singleStr[[i]]$set_seqSt_update_neighbSt(j, private$singleStr[[i]]$get_seqSt_leftneighb(index = j))
           if (testing){
             event_acceptance[k] <- TRUE
             r_jk[k] <- r
             applied_event[k] <- paste("SSEc", "left", sep = "_")
           } 
         } else {
           ## copy from right neighbor
           private$singleStr[[i]]$set_seqSt_update_neighbSt(j, private$singleStr[[i]]$get_seqSt_rightneighb(index = j))
           if (testing){
             event_acceptance[k] <- TRUE
             r_jk[k] <- r
             applied_event[k] <- paste("SSEc", "right", sep = "_")
           } 
         }
      } else if (testing) {
        r_jk[k] <- r
      }
    }
  }
  if(testing){
    list(CFTP_chosen_singleStr = private$CFTP_chosen_singleStr,
         CFTP_chosen_site = private$CFTP_chosen_site,
         CFTP_event = private$CFTP_event,
         CFTP_random = private$CFTP_random,
         event_acceptance = event_acceptance,
         applied_event = applied_event,
         r_jk = r_jk,
         r_m = private$CFTP_highest_rate)
  }
})

## @description
## Public Method. Applies the CFTP algorithm.
## 
## @param steps minimum number of steps to apply (default 10000)
## @param testing default FALSE. TRUE for testing output
## 
## @return combiStructureGenerator instance when testing FALSE. Testing output when testing TRUE.
combiStructureGenerator$set("public", "cftp", function(steps = 10000, testing = FALSE) {
    # Set a variable to track when the $seq of the 2 combi instances become equal
    equal <- FALSE
    while(!equal) {
      # Sample the CFTP steps 
      self$cftp_event_generator(steps)
      # Copy (deep clone) the combiStructure instance to generate 2 instances
      combi_u <- self$copy()
      combi_m <- self$copy()
      for(str in 1:length(private$singleStr)){
        # Set the sequences for each as all m states and all u states
        combi_u$get_singleStr(str)$cftp_all_equal(state = "U")
        combi_m$get_singleStr(str)$cftp_all_equal(state = "M")
      }
       for(str in 1:length(private$singleStr)){
        # Update the neighbSt according to the new sequence
        combi_u$get_singleStr(str)$init_neighbSt()
        combi_m$get_singleStr(str)$init_neighbSt()
        ### note that rate trees are not updated in combi_u and combi_m 
      }
        combi_u$cftp_apply_events() 
        combi_m$cftp_apply_events() 
      
        equal_str <- c()
        for(str in 1:length(private$singleStr)){
          equal_str[str] <- all(combi_u$get_singleStr(str)$get_seq() == combi_m$get_singleStr(str)$get_seq())
        }
        equal <- all(equal_str)
        steps <- 2*steps
    }
    for(str in 1:length(private$singleStr)){
      # update converged combiStructure ratetree
      combi_u$get_singleStr(str)$initialize_ratetree()
    }
    if(testing){
      return(list(combi_u = combi_u,
                  combi_m = combi_m,
                  total_steps = length(private$CFTP_event)))
    } else {
      return(combi_u)
    }
    
})






