
############ Validation of Vectors containing frequency values #################

test_that("validate_freqVectorSums1", {
  # Test valid input
  valid_freqVector <- c(0.2, 0.3, 0.5)
  expect_true(validate_freqVectorSums1(valid_freqVector)$valid)

  # Test invalid input - non-numeric freqVector
  non_numeric_freqVector <- c(0.2, "invalid", 0.5)
  expect_error(validate_freqVectorSums1(non_numeric_freqVector),
               info = "funtion does not return an error for non-numeric vector input")
  non_vector_freqVector <- matrix(c(.2, .3, .5), nrow = 1)
  expect_error(validate_freqVectorSums1(non_vector_freqVector),
               info = "function does not return an error for non-vector input")

  # Test function reports non-valid state when given freqs don't sum to 1
  freqVector_notSum1 <- c(0.3, 0.3, 0.5)
  expect_false(validate_freqVectorSums1(freqVector_notSum1)$valid,
               info = "fails reporting non-valid freqs state when freqs dont sum 1")

  # Function reports valid state when given freqs sum to 1
  expect_true(validate_freqVectorSums1(valid_freqVector)$valid,
              info = "fails reporting valid freqs state when freqs sum 1")

})

test_that("validate_freqVector_elements", {
  # Test valid input
  valid_freqVector <- c(0.2, 0.3, 0.5)
  expect_no_error(validate_freqVector_elements(valid_freqVector))

  # Test invalid input - non-numeric freqVector
  non_numeric_freqVector <- c(0.2, "invalid", 0.5)
  expect_error(validate_freqVector_elements(non_numeric_freqVector),
               info = "funtion does not return an error for non-numeric vector input")
  non_vector_freqVector <- matrix(c(.2, .3, .5), nrow = 1)
  expect_error(validate_freqVector_elements(non_vector_freqVector),
               info = "function does not return an error for non-vector input")

  # Test function reports non-valid state when given freqs elements > 1 or < 0
  freqVector_invalid <- c(-0.3, 0.3, 0.5)
  expect_false(validate_freqVector_elements(freqVector_invalid)$valid,
               info = "fails reporting non-valid freqs state when elements arent between 0 and 1")

  # Function reports valid state when given freqs sum to 1
  valid_freqVector <- c(0, 0, 1)
  expect_true(validate_freqVector_elements(valid_freqVector)$valid,
              info = "fails reporting valid freqs state when elements are between 0 and 1")

})

test_that("listFreqVector_validation", {
  invalid_vector_list <- list(c(1, -1, 1, 1), "a")
  valid_vector_list <- list(c(1, -1, 1, 1), c(1, 2, 3))
  valid_listName <- "Qci"
  wrong_listName <- 1
  # Function throws an error when element of list is not a matrix
  expect_error(listFreqVector_validation(valid_vector_list, wrong_listName),
               info = "fails throwing error for non character string listName argument")
  # Function returns correct listName
  expect_equal(listFreqVector_validation(valid_vector_list, valid_listName)$listName, valid_listName,
               info = "listName output is not equal to given listName as argument")
  # Function throws an error when element of list is not a matrix
  expect_error(listFreqVector_validation(invalid_vector_list, valid_listName),
               info = "fails throwing error for non-vector list listTransitionMatrix argument")
  # Function returns list of validationStates with number of elements to input list
  expect_true(length(listFreqVector_validation(valid_vector_list, valid_listName)$validationStates) == 2,
              info = "validationStates length different form provided listFreqVector length")
})

#########################################

explore_island_n <- c(5,10)
explore_mu_values <- c(0.1, 1, 2)
reps_for_T1error <- 1000
branchLength <- 1
rep_n <- 100 
for (island in 1:length(explore_island_n)){
  island_n <- explore_island_n[island]
  infoStr <- data.frame(n = rep(10, island_n),
                        globalState = rep("U", island_n))
  for (mu in 1:length(explore_mu_values)){
    custom_params <- get_parameterValues()
    custom_mu <- explore_mu_values[mu]
    custom_params$mu <- custom_mu
    print(paste("######### Explored case: mu =", custom_mu, ", island_n =", island_n))
    obj <- combiStructureGenerator$new(infoStr, params = custom_params)
    rejecH0_IWEnumber_n <- 0
    rejecH0_IWEisland_n <- 0
    for (i in 1:reps_for_T1error){
      number_IWEs <- c()
      island_IWEs <- list()
      for (rep in 1:rep_n){
        output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
        # I expect 10 IWE events (10 islands * rate/island and bL unit = 1)
        if (output$IWE_event){
          expect_equal(length(output$IWE_times), length(output$islands), 
                       info = paste("Different number of islands sampled for the number of IWEs sampled in rep_n: ", rep))
          number_IWEs[rep] <- length(output$IWE_times)
          island_IWEs[[rep]] <- output$islands
        } else {
          # Expected number of IWEs is 10. Quite unlikely to sample none
          number_IWEs[rep] <- 0
          island_IWEs[[rep]] <- NA
        }
      }
      # T-test for expected number of IWEs
      if(t.test(number_IWEs - custom_mu*island_n)$p.value <= 0.05){
        rejecH0_IWEnumber_n <- rejecH0_IWEnumber_n +1
      }
      
      # Chi-squared test for uniformity of sampled islands
      sampled_islands <- unlist(island_IWEs)
      sampled_islands <- sampled_islands[!is.na(sampled_islands)]
      # Ensure all indices from 1 to 10 are included
      sampled_islands <- factor(sampled_islands, levels = 1:island_n)
      expected_freq <- rep(length(sampled_islands)/island_n, island_n)
      observed_freq <- table(sampled_islands)
      if (unique(expected_freq)>=5){
        chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n))
      } else {
        chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n), simulate.p.value = TRUE, B = 10000)
      }
      if (chi_square_test$p.value <= 0.05){
        rejecH0_IWEisland_n <- rejecH0_IWEisland_n +1 
      }
    }
    print(paste("Proportion of times Observed number of IWEs significantly deviates from expected:", rejecH0_IWEnumber_n/reps_for_T1error))
    print(paste("Proportion of times IWEs are assigned to islands in a significantly different way from uniform:", rejecH0_IWEisland_n/reps_for_T1error))
  }
}


## Missing to test for the 1 island case:
#####################################################################################################################
island_n <- 10
infoStr <- data.frame(n = rep(10, island_n),
                      globalState = rep("U", island_n))
custom_params <- get_parameterValues()
custom_mu <- 0.01
custom_params$mu <- custom_mu
obj <- combiStructureGenerator$new(infoStr, params = custom_params)
branchLength <- 1
rep_n <- 100
number_IWEs <- c()
island_IWEs <- list()
for (rep in 1:rep_n){
  output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
  # I expect 10 IWE events (10 islands * rate/island and bL unit = 1)
  if (output$IWE_event){
    expect_equal(length(output$IWE_times), length(output$islands), 
                 info = paste("Different number of islands sampled for the number of IWEs sampled in rep_n: ", rep))
    number_IWEs[rep] <- length(output$IWE_times)
    island_IWEs[[rep]] <- output$islands
  } else {
    # Expected number of IWEs is 10. Quite unlikely to sample none
    number_IWEs[rep] <- 0
    island_IWEs[[rep]] <- NA
  }
}
mean(number_IWEs)
table(unlist(island_IWEs))

# T-test for expected number of IWEs
t_test <- t.test(number_IWEs - custom_mu*island_n)

# Chi-squared test for uniformity of sampled islands
sampled_islands <- unlist(island_IWEs)
sampled_islands <- sampled_islands[!is.na(sampled_islands)]
# Ensure all indices from 1 to 10 are included
sampled_islands <- factor(sampled_islands, levels = 1:island_n)
expected_freq <- rep(length(sampled_islands)/island_n, island_n)
observed_freq <- table(sampled_islands)
if (unique(expected_freq)>=5){
  chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n))
} else {
  chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n), simulate.p.value = TRUE, B = 10000)
}


############################

