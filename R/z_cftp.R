singleStructureGenerator$set("public", "get_Q", function() {
  private$Q
})

combiStructureGenerator$set("public", "cftp", function() {
  
  highest_rate <- 0
  singleStr_n <- c()
  for(str in 1:length(private$singleStr)){
    highest_rate <- max(highest_rate, 
                        max(unlist(private$singleStr[[str]]$get_Q())))
    singleStr_n[str] <- length(private$singleStr[[str]]$get_seq())
  }
  
  # Propose 1 site and what may happen to it
  chosen_singleStr <- sample(1:length(private$singleStr), 1, prob=singleStr_n)
  chosen_site <- sample(1:length(private$singleStr[[str]]$get_seq()), 1)
  print(chosen_singleStr)
  print(chosen_site)
  event <- sample(1:3, 1)  ## 1, 2, 3: go to u, p, m
  
  random_threshold <- runif(1)
  
    ## next steps and todos:
    ##
    ## for n steps (corresponding to time steps from -n, -n+1, ... -1; e.g for an initial n=10000)
    ## sample a chosen site, an event and a random threshold; store them either in vectors
    ## (where time step -k is stored at index k) or maybe even in a specific R6Class object.
    ## (Or turn it into and R6Object later as part of refactoring.)
    ## To be decided: are the two CpG sequences represented by two separate combiStructures or
    ## by two simpler representations of {u,p,m}-Sequences, e.g. vectors.
    ##
    ## Apply these steps, depending on probabilities, beginning at time -n to two
    ## CpG sequences, one consisting of all sites u and the other with all sites m.
    ## An event is applied iff it probability is greather than the random_threshold
    ## (by the way, iff is math slang for "if and only if")
    ##
    ## If at time 0, that is after all n possible events, the two CpG sequences agree in each position
    ## this sequence is the result. Otherwise extend the series of events by simulating chosen sites,
    ## events and a random thresholds for time points -2n, -2n+1,... -n-1. (Important: the events etc.
    ## for -n, -n+1, ...,-1 are still the ones from above.) Then continue applying the steps as above
    ## with n replaced by 2n to two sequences u,u,u,...,u and m,m,m,...,m starting at time -2n. 
    ##
    ## Continue as above until with always doubling n until the two sequences end in the same state
    ## at time 0.
    ##
    ## Both of us should check whether this is really the CFTP algorithm and understand why the CFTP algorithm
    ## leads to a sequence sampled from the equilibrium distribution of the Castillo-Vicente--Metzler model.
    ##
    ## Testing: Unit tests should check whether applications of events as above to two sequences a and b
    ## of which a <= b is fulfilled ("a <= b" means that "a < b" or "a = b"), alwas lead to sequences a' and b'
    ## with a' <= b'.
    ##
    ## More time-consuming function tests could check with certain summary statistics, e.g. combinations of
    ## neighboring pairs or triples within CpGs, in non-CpG or at boundaries have the same frequencies in initial
    ## cftp-generated states as in states that evolved for a while after stating in such a state.
    ##
    ## Extend CFTP by allowing also IWEs.
    ##
    ## Perhaps later: look for possibilities to improve structure of implementation (refactoring) and efficiency
   
})



