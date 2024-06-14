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
})



