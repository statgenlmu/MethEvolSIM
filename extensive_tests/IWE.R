## This test checks whether the number of IWEs simulated is in a plausible range for a given rate of IWEs (mu)
## and whether the IWEs are assigned to the islands with equal probability
## To run it the following commands need to be used:
library(devtools)
load_all()

## Test for the Expected number of IWEs and the distribution of IWEs across islands

# Explore cases with the following number of islands:
explore_island_n <- c(5,10)
# Explore cases with the following IWE rates (mu):
explore_mu_values <- c(0.1, 1, 2)
# Set the number of repetitions to calculate the proportion of Type 1 error (significant deviations from expected):
reps_for_T1error <- 1000
# Set branch length for simulations
branchLength <- 1
# Set the number of simulations conducted per combination of number of islands and IWE rate
rep_n <- 100 
for (island in 1:length(explore_island_n)){
  # Generate the structure containing the corresponding number of islands
  island_n <- explore_island_n[island]
  infoStr <- data.frame(n = rep(10, island_n),
                        globalState = rep("U", island_n))
  for (mu in 1:length(explore_mu_values)){
    # Set the corresponding mu value
    custom_params <- get_parameterValues()
    custom_mu <- explore_mu_values[mu]
    custom_params$mu <- custom_mu
    print(paste("######### Explored case: mu =", custom_mu, ", island_n =", island_n))
    # Generate the combiStructure
    obj <- combiStructureGenerator$new(infoStr, params = custom_params)
    # Initialize a counter for the number of times the observd value significantly deviates from expected
    rejecH0_IWEnumber_n <- 0
    rejecH0_IWEisland_n <- 0
    for (i in 1:reps_for_T1error){
      # Initialize the objects to store the number of IWEs and islands affected 
      number_IWEs <- c()
      island_IWEs <- list()
      for (rep in 1:rep_n){
        output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
        if (output$IWE_event){ # If at least 1 IWE was sampled:
          if (length(output$IWE_times)!= length(output$islands)){
            print(paste("Different number of islands sampled for the number of IWEs sampled in rep_n: ", rep))
            print("IWE_times:")
            print(output$IWE_times)
            print("islands:")
            print(output$islands)
          }
          number_IWEs[rep] <- length(output$IWE_times)
          island_IWEs[[rep]] <- output$islands
        } else { # If no IWE was sampled
          number_IWEs[rep] <- 0
          island_IWEs[[rep]] <- NA
        }
      }
      # Poisson test for expected number of IWEs
      if(poisson.test(sum(number_IWEs), rep_n*custom_mu*island_n)$p.value <= 0.05){
        rejecH0_IWEnumber_n <- rejecH0_IWEnumber_n +1
      }
      
      # Chi-squared test for uniformity of sampled islands
      sampled_islands <- unlist(island_IWEs)
      sampled_islands <- sampled_islands[!is.na(sampled_islands)]
      # Ensure all indices from 1 to 10 are included
      sampled_islands <- factor(sampled_islands, levels = 1:island_n)
      expected_freq <- rep(length(sampled_islands)/island_n, island_n)
      observed_freq <- table(sampled_islands)
      if (unique(expected_freq)>=5){ # If assumptions for Chi-squared are met:
        chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n))
      } else { # If minimum frequencies are not met: 
        chi_square_test <- chisq.test(observed_freq, p = rep(1/island_n, island_n), simulate.p.value = TRUE, B = 10000)
      }
      if (chi_square_test$p.value <= 0.05){
        rejecH0_IWEisland_n <- rejecH0_IWEisland_n +1 
      }
    }
    # Print the result for the combination of mu and island number
    print(paste("Proportion of times Observed number of IWEs significantly deviates from expected:", rejecH0_IWEnumber_n/reps_for_T1error))
    print(paste("Proportion of times IWEs are assigned to islands in a significantly different way from uniform:", rejecH0_IWEisland_n/reps_for_T1error))
  }
}