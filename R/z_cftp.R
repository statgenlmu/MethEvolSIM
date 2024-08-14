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

singleStructureGenerator$set("public", "get_leftneighbSt", function(i){
  floor((self$get_neighbSt(i)-1) / 3) +1
})

singleStructureGenerator$set("public", "get_rightneighbSt", function(i){
  (self$get_neighbSt(i)-1) %% 3 +1
})

## Make a singleStructure with the same segment lengths and parameters
## as the focal one but all states are m
singleStructureGenerator$set("public", "cftp_all_equal", function(state, testing = FALSE) {
  n <- length(private$seq)
  if (state == "U"){
    private$globalState <- "U"
    private$seq <- rep(1L, n)
  }
  if (state == "M"){
    private$globalState <- "M"
    private$seq <- rep(3L, n)
  }
})

##singleStructureGenerator$set("public", "get_siteR", function(index = NULL){
##  if(is.null(index)){
##    private$siteR
##  } else{
##    private$siteR[index]
##  }
##})

##singleStructureGenerator$set("public", "get_neighbSt", function(index = NULL){
##  if(is.null(index)){
##    private$neighbSt
##  } else{
##    private$neighbSt[index]
##  }
##})


## set sequence, update neighbors BUT NOT THE ratetree
singleStructureGenerator$set("public", "set_seq", function(index, newSt){
  private$seq[index] <- newSt
  private$update_neighbSt(index)
})


## note that this does not update the ratetrees; needs to be done after the whole CFTP run
combiStructureGenerator$set("public", "cftp_apply_events", function() {
  
  for(k in length(private$CFTP_event):1) {
    ### applies the CFTP_events from -n to 1 to the combiStructure
    i <- private$CFTP_chosen_singleStr[k]
    j <- private$CFTP_chosen_site[k]
    siteR <- private$singleStr[[i]]$get_siteR(j)
    neighbSt <- private$singleStr[[i]]$get_neighbSt(j)
    oldSt <- private$singleStr[[i]]$get_seq()[j]
    if(private$CFTP_event[k] < 4) {
      newSt <- private$CFTP_event[k]      
      if(oldSt != newSt) {
        r <- private$singleStr[[i]]$get_Qi(siteR = siteR, oldSt = oldSt, newSt = newSt)
        if( r/private$CFTP_highest_rate > private$CFTP_random[k] ) {
            private$singleStr[[i]]$set_seq(j, newSt)
        }
      }
    } else {
      r <- (1-private$singleStr[[i]]$get_iota())/2
      if( r/private$CFTP_highest_rate > private$CFTP_random[k] ) {
         if(private$CFTP_event[k] == 4) {
           ##copy from left neighbor
           private$singleStr[[i]]$set_seq(j, private$singleStr[[i]]$get_leftneighbSt())
         } else {
           ## copy from right neighbor
           private$singleStr[[i]]$set_seq(j, private$singleStr[[i]]$get_leftneighbSt())
         }
      }
    }
  }
})

combiStructureGenerator$set("public", "cftp", function() {
    equal <- FALSE
    steps <- 10000
    while(!equal) {
        self$cftp_event_generator(steps)
        combi_u <- self$clone()
        combi_m <- self$clone()
        for(str in 1:length(private$singleStr)){
            combi_u$get_singleStr(str)$cftp_all_equal(state = "U")
            combi_m$get_singleStr(str)$cftp_all_equal(state = "M")
            ### note that rate trees are not updated in combi_u and combi_m
            combi_u$get_singleStr(str)$init_neighbSt()
            combi_m$get_singleStr(str)$init_neighbSt()
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
      combi_u$get_singleStr(str)$initialize_ratetree()
    }
    return(combi_u)
})


combiStructureGenerator$set("public", "cftp_event_generator", function(steps) {
  
  private$CFTP_highest_rate <- 0
  singleStr_n <- c()
  for(str in 1:length(private$singleStr)){
    private$CFTP_highest_rate <- max(private$CFTP_highest_rate,            
                         max(c(unlist(private$singleStr[[str]]$get_Qi()), 
                               (1-private$singleStr[[str]]$get_iota())/2)))
    singleStr_n[str] <- length(private$singleStr[[str]]$get_seq())
  }
  chosen_singleStr <- integer(length=steps)
  chosen_site <- integer(length=steps)
  event <- integer(length=steps)
  random_threshold <- numeric(length=steps)
  for(n in steps:1) {
    # For generation -n
    # Propose 1 site and what may happen to it
    chosen_singleStr[n] <- sample(1:length(private$singleStr), 1, prob=singleStr_n)
    chosen_site[n] <- sample(1:length(private$singleStr[[str]]$get_seq()), 1)
    event[n] <- sample(1:5, 1)  ## 1,..5: go to u, p, m by SSEi or copy left, copy right.                        
  
    random_threshold[n] <- runif(1)
  }
  private$CFTP_chosen_singleStr <- c(private$CFTP_chosen_singleStr, chosen_singleStr)
  private$CFTP_chosen_site <- c(private$CFTP_chosen_site, chosen_site)
  private$CFTP_event <- c(private$CFTP_event, event)
  private$CFTP_random <- c(private$CFTP_random, random_threshold) 
})



