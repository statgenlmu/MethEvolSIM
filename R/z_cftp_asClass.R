combiStructureGenerator$set("public", "get_highest_rate", function() {
  CFTP_highest_rate <- 0 # ensure minimum value of 0
  for(str in 1:length(private$singleStr)){
    CFTP_highest_rate <- max(CFTP_highest_rate,            
                                     max(c(unlist(private$singleStr[[str]]$get_Qi()), 
                                           (1-private$singleStr[[str]]$get_iota())/2)))
    
  }
  CFTP_highest_rate
})

combiStructureGenerator$set("public", "get_singleStr_siteNumber", function() {
  singleStr_site_n <- c()
  for(str in 1:length(private$singleStr)){
    singleStr_site_n[str] <- length(private$singleStr[[str]]$get_seq())
    
  }
  singleStr_site_n
})

# Method self$get_singleStr_number() exists

combiStructureGenerator$set("public", "get_CFTP_info", function() {
  private$CFTP_info
})

combiStructureGenerator$set("public", "set_CFTP_info", function(CFTP_instance) {
  if(class(CFTP_instance)[1] != "cftpStepGenerator") stop("Argument 'CFTP_instance' must be of class cftpStepGenerator")
  private$CFTP_info <- CFTP_instance
})

cftpStepGenerator <- R6::R6Class("cftpStepGenerator",
                                 public = list(
                                   ##TODO: document combi instance attributes
                                   singleStr_number = NULL,
                                   singleStr_siteNumber = NULL,
                                   ##TODO: document CFTP attributes
                                   CFTP_highest_rate = 0,
                                   CFTP_chosen_singleStr = integer(length=0),
                                   CFTP_chosen_site = integer(length=0),
                                   CFTP_event = integer(length=0),
                                   CFTP_random = numeric(length=0),
                                   
                                   ##TODO: document
                                   generate_events = function(steps = 10000, testing = FALSE) {
                                     if(!(is.numeric(steps) && length(steps) == 1 && steps == floor(steps) && steps >= 1)){
                                       stop("'index' must be one number >=1 without decimals")
                                     }
                                     
                                     # Generate for each CFTP step the event and location to apply it and a threshold for acceptance/rejection
                                     old_steps <- length(self$CFTP_event) # Check the number of already existing steps
                                     if (steps <= old_steps){stop("The given number of steps has already been generated")}
                                     
                                     # Get the number of new steps
                                     new_steps <- steps - old_steps 
                                     
                                     # Initialize the vectors to store the CFTP info for the new steps
                                     chosen_singleStr <- integer(length=new_steps)
                                     chosen_site <- integer(length=new_steps)
                                     event <- integer(length=new_steps)
                                     random_threshold <- numeric(length=new_steps)
                                     
                                     # For each CFTP step (from past to present)
                                     for(n in new_steps:1) {
                                       # For generation -n
                                       # Propose 1 site and what may happen to it
                                       chosen_singleStr[n] <- sample(1:self$singleStr_number, 1, prob=self$singleStr_siteNumber)
                                       chosen_site[n] <- sample(1:self$singleStr_siteNumber[chosen_singleStr[n]], 1)
                                       event[n] <- sample(1:5, 1)  ## 1,2,3: go to u, p, m by SSEi ## 4,5: copy left, copy right.                        
                                       # Sample a threshold to accept or reject event
                                       random_threshold[n] <- runif(1) # numerical value between 0 and 1
                                     }
                                     
                                     # Add the new step info to the existing ones
                                     self$CFTP_chosen_singleStr <- c(self$CFTP_chosen_singleStr, chosen_singleStr)
                                     self$CFTP_chosen_site <- c(self$CFTP_chosen_site, chosen_site)
                                     self$CFTP_event <- c(self$CFTP_event, event)
                                     self$CFTP_random <- c(self$CFTP_random, random_threshold)
                                     
                                     if(testing){
                                       list(CFTP_chosen_singleStr = self$CFTP_chosen_singleStr,
                                            CFTP_chosen_site = self$CFTP_chosen_site,
                                            CFTP_event = self$CFTP_event,
                                            CFTP_random = self$CFTP_random)
                                     }
                                     
                                   },
                                   
                                   ##TODO: document 
                                   initialize = function(singleStr_number, singleStr_siteNumber, CFTP_highest_rate) {
                                     # Initialize a new instance with the info of the corresponding combiStrucutre instance
                                     self$singleStr_number <- singleStr_number
                                     self$singleStr_siteNumber <- singleStr_siteNumber
                                     # Set CFTP_highest_rate to be the highest rate across all singleStr withing combiStr instance
                                     self$CFTP_highest_rate <- CFTP_highest_rate
                                     
                                   }
                                 ))