# TODO: Ensure the number of steps in a single vector remains
# below the number of elements in a numeric vector of size 1.2 GB.#

cftpStepGenerator_new <- R6::R6Class("cftpStepGenerator_new",
                                     public = list(
                                       #' @field singleStr_number Public attribute: Number of singleStr instances
                                       singleStr_number = NULL,
                                       #' @field singleStr_siteNumber Public attribute: Number of sites in singleStr instances
                                       singleStr_siteNumber = NULL,
                                       #' @field CFTP_highest_rate Public attribute: CFTP highest rate
                                       CFTP_highest_rate = 0,
                                       
                                       #' @description
                                       #' Create a new instance of class cftpStepGenerator with the info of the corresponding combiStrucutre instance
                                       #'
                                       #' @param singleStr_number Number of singleStr instances
                                       #' @param singleStr_siteNumber Number of sites in singleStr instances
                                       #' @param CFTP_highest_rate CFTP highest rate across all singleStr withing combiStr instance 
                                       #' 
                                       #' @return A new instance of class cftpStepGenerator
                                       initialize = function(singleStr_number, singleStr_siteNumber, CFTP_highest_rate) {
                                         # Initialize a new instance with the info of the corresponding combiStrucutre instance
                                         self$singleStr_number <- singleStr_number
                                         self$singleStr_siteNumber <- singleStr_siteNumber
                                         # Set CFTP_highest_rate to be the highest rate across all singleStr withing combiStr instance
                                         self$CFTP_highest_rate <- CFTP_highest_rate
                                         
                                       },
                                       
                                       ##TODO: document
                                       number_steps = 0,
                                       ##TODO: update documentation for lists
                                       #' @field CFTP_chosen_singleStr Public attribute: chosen singleStr index at each CFTP step
                                       CFTP_chosen_singleStr = list(),
                                       #' @field CFTP_chosen_site Public attribute: chosen site index at each CFTP step
                                       CFTP_chosen_site = list(),
                                       #' @field CFTP_event Public attribute: type of CFTP event at each CFTP step.
                                       #' @description
                                       #' 1: SSEi to unmethylated, 2: SSEi to partially methylated, 3: SSEi to methylated
                                       #' 4: SSEc copy left state, 5: SSEc copy right state
                                       CFTP_event = list(),
                                       #' @field CFTP_random Public attribute: CFTP threshold at each CFTP step
                                       CFTP_random = list(),
                                       ##TODO: document
                                       steps_perVector = NULL,
                                       
                                       ##TODO: Update documentation
                                       #' @description
                                       #' Public Method. Generates the events to apply for CFTP.
                                       #' 
                                       #' @param steps Integer value >=1
                                       #' @param testing default FALSE. TRUE for testing output
                                       #' 
                                       #' @return NULL when testing FALSE. Testing output when testing TRUE.
                                       #'
                                       #' @details
                                       #' The function add steps to the existing ones. 
                                       #' If called several times the given steps need to be higher than the sum of steps generated before.
                                       generate_events = function(steps = 10000, testing = FALSE){
                                         if(!(is.numeric(steps) && length(steps) == 1 && steps == floor(steps) && steps >= 1)){
                                           stop("'index' must be one number >=1 without decimals")
                                         }
                                         
                                         # Set the number of steps per vector
                                         if(is.null(self$steps_perVector)){
                                           self$steps_perVector <- steps
                                         }
                                         
                                         # Check the number of already existing steps
                                         old_steps <- self$number_steps
                                         if (steps <= old_steps){stop("The given number of steps has already been generated")}
                                         
                                         # Get the number of new steps
                                         new_steps <- steps - old_steps
                                         
                                         # Add the info of the new steps to the existing ones
                                         # handling potential errors (e.g. memmory limits)
                                         tryCatch({
                
                                           # Set a variable to store the number of remaining steps to sample 
                                           remaining_steps <- new_steps
                                           
                                           while(remaining_steps > 0){
                                             
                                             # Initialize the vectors to store the CFTP info for the new steps
                                             chosen_singleStr <- integer(length=self$steps_perVector)
                                             chosen_site <- integer(length=self$steps_perVector)
                                             event <- integer(length=self$steps_perVector)
                                             random_threshold <- numeric(length=self$steps_perVector)
                                             
                                             # Generate for each CFTP step the event and location to apply it and a threshold for acceptance/rejection
                                             # For each CFTP step (from past to present)
                                             for(n in self$steps_perVector:1) {
                                               # For generation -n
                                               # Propose 1 site and what may happen to it
                                               chosen_singleStr[n] <- sample(1:self$singleStr_number, 1, prob=self$singleStr_siteNumber)
                                               chosen_site[n] <- sample(1:self$singleStr_siteNumber[chosen_singleStr[n]], 1)
                                               event[n] <- sample(1:5, 1)  ## 1,2,3: go to u, p, m by SSEi ## 4,5: copy left, copy right.                        
                                               # Sample a threshold to accept or reject event
                                               random_threshold[n] <- runif(1) # numerical value between 0 and 1
                                             }
                                             
                                             # Set the current list index to store CFTP info 
                                             current_listIndex <- length(self$chosen_singleStr) + 1
                                             
                                             # Add the info of the new steps to the existing ones
                                             self$CFTP_chosen_singleStr[[current_listIndex]] <- c(self$CFTP_chosen_singleStr, chosen_singleStr)
                                             self$CFTP_chosen_site[[current_listIndex]] <- c(self$CFTP_chosen_site, chosen_site)
                                             self$CFTP_event[[current_listIndex]] <- c(self$CFTP_event, event)
                                             self$CFTP_random[[current_listIndex]] <- c(self$CFTP_random, random_threshold)
                                             
                                             
                                             # Update the number of remaining steps to sample 
                                             remaining_steps <- remaining_steps - self$steps_perVector
                                           }
                                           
                                         }, error = function(e) {
                                           message("Error encountered: ", conditionMessage(e))
                                           
                                           
                                           warning(paste0("Reducing the number of steps to", steps, "the number of steps per vector to", self$steps_perVector, "and retrying."))
                                           
                                           # If there is at least one more step retry, if not return something that can be catched in the outer function
                                           if(new_steps < 1 || self$steps_perVector < 1) {
                                             warning("Failure to increase the number of CFTP steps")
                                             return("failed")
                                           } else {
                                             # Retry
                                             return(self$generate_events(steps = old_steps + new_steps, testing))
                                           }
                                         })

                                         
                                         # Update the number of already existing steps
                                         self$number_steps <- steps
                                         

                                         if(testing){
                                           list(old_steps = old_steps,
                                                new_steps = new_steps,
                                                current_listIndex = current_listIndex)
                                         }
                                       }
                                     ))