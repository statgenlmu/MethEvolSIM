test_that("get_MeanFreqP structures with equal length islands and non-islands", {
  # Two tips and two islands
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(0, 0, 0), c(0, 0, 0)), # tip 1
    list(c(0, 0, 0), c(1, 1, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 0,
               info = "incorrect mean in two tips and two islands, freq 0")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 0,
               info = "incorrect mean in two tips and two non-islands, freq 0")
  
  # Two tips and two islands
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(.5, .5, 0), c(0, 0, .5)), # tip 1
    list(c(0, .5, .5), c(1, .5, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), .5,
               info = "incorrect mean in two tips and two islands, freq .5")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), .5,
               info = "incorrect mean in two tips and two non-islands, freq .5")
  
  # Two tips and 1 island
  index_islands <- c(2) # island index 2
  index_nonislands <- c(2) # non-island index 2
  data <- list(
    list(c(.5, .5, 0), c(0, 0, .5)), # tip 1
    list(c(0, .5, .5), c(1, .5, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 1/3,
               info = "incorrect mean in two tips and one islands, freq 1/3")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 1/3,
               info = "incorrect mean in two tips and one non-island, freq 1/3")
  index_islands <- c(1) # island index 1
  index_nonislands <- c(1) # island index 1
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 2/3,
               info = "incorrect mean in two tips and one islands, freq 2/3")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 2/3,
               info = "incorrect mean in two tips and one non-island, freq 2/3")
  
  # One tip, one island
  index_islands <- c(1)
  index_nonislands <- c(1)
  sample_n <- 1
  
  data <- list(c(0.5, 1, 0.5)) # single tip
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 2/3,
               info = "incorrect mean in one tip one island")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 2/3,
               info = "incorrect mean in one tip one non-island")

  data <- list(c(1, 1, 0)) # single tip, no 0.5 state
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 0,
               info = "incorrect mean in one tip one island no .5 state")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 0,
               info = "incorrect mean in one tip one non-island no .5 state")
})

test_that("get_MeanFreqP structures with different length (islands and non-islands)", {
  # Two tips and two islands, freq 0
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(0, 0, 0, 0, 0, 0, 0), c(0, 0, 0)), # tip 1
    list(c(0, 0, 0, 0, 0, 0, 0), c(1, 1, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), 0,
               info = "incorrect mean in two tips and two islands, freq 0")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), 0,
               info = "incorrect mean in two tips and two non-islands, freq 0")
  
  # Two tips and two islands, freq .5 evenly distributed
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(.5, .5, 0, 0, 0, .5), c(.5, 0, 0, .5)), # tip 1
    list(c(0, .5, .5, 1, 1, .5), c(1, .5, 1, .5)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), .5,
               info = "incorrect mean in two tips and two islands, freq .5")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), .5,
               info = "incorrect mean in two tips and two non-islands, freq .5")
  
  # Two tips and two islands, freq .5 unevenly distributed (fist island mean 3/4, second island mean 1/4)
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(.5, .5, 0, .5, .5, .5), c(.5, 0, 0, 0)), # tip 1
    list(c(0, .5, .5, 1, .5, .5), c(1, 0, 1, .5)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqP(index_islands, data, sample_n), .5,
               info = "incorrect mean in two tips and two islands, freq .5")
  expect_equal(get_nonislandMeanFreqP(index_nonislands, data, sample_n), .5,
               info = "incorrect mean in two tips and two non-islands, freq .5")
 
})


test_that("get_MeanFreqM structures with equal length islands and non-islands", {
  # Two tips and two islands
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(0, 0, 0), c(0, 0, 0)), # tip 1
    list(c(0, 0, 0), c(.5, .5, .5)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 0,
               info = "incorrect mean in two tips and two islands, freq 0")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 0,
               info = "incorrect mean in two tips and two non-islands, freq 0")
  
  # Two tips and two islands
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(1, 1, 0), c(0, 0, .5)), # tip 1
    list(c(1, 1, .5), c(1, .5, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 0.5,
               info = "incorrect mean in two tips and two islands, freq 0.5")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 0.5,
               info = "incorrect mean in two tips and two non-islands, freq 0.5")
  
  # Two tips and 1 island
  index_islands <- c(2) # island index 2
  index_nonislands <- c(2) # non-island index 2
  data <- list(
    list(c(.5, .5, 0), c(0, 0, .5)), # tip 1
    list(c(0, .5, 1), c(1, 1, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), .5,
               info = "incorrect mean in two tips and one island, freq 0.5")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), .5,
               info = "incorrect mean in two tips and one non-island, freq 0.5")
  
  # One tip, one island
  index_islands <- c(1)
  index_nonislands <- c(1)
  sample_n <- 1
  
  data <- list(c(1, 1, 0.5)) # single tip
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 2/3,
               info = "incorrect mean in one tip one island")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 2/3,
               info = "incorrect mean in one tip one non-island")
  
  data <- list(c(0, 0, 0)) # single tip, no methylated state
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 0,
               info = "incorrect mean in one tip one island no methylated state")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 0,
               info = "incorrect mean in one tip one non-island no methylated state")
})

test_that("get_MeanFreqM structures with different length (islands and non-islands)", {
  # Two tips and two islands, freq 1/4
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(0, 0, 0, 0, 0, 0, 0), c(0, 0, 0)), # tip 1
    list(c(0, 0, 0, 0, 0, 0, 0), c(1, 1, 1)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 1/4,
               info = "incorrect mean in two tips and two islands, freq 1/4")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 1/4,
               info = "incorrect mean in two tips and two non-islands, freq 1/4")
  
  # Two tips and two islands, freq .5 evenly distributed
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(.5, .5, 1, 1, 1, .5), c(1, 1, .5, 0)), # tip 1
    list(c(0, 1, .5, 1, 1, .5), c(1, 1, 0, 0)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), 0.5,
               info = "incorrect mean in two tips and two islands, freq .5")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), 0.5,
               info = "incorrect mean in two tips and two non-islands, freq .5")
  
  # Two tips and two islands, freq .5 unevenly distributed (first island mean 3/4, second island mean 1/4)
  index_islands <- c(1, 2)
  index_nonislands <- c(1, 2)
  data <- list(
    list(c(1, .5, 0, 1, 1, 1), c(.5, 1, 0, .5)), # tip 1
    list(c(1, 1, .5, 1, 1, 1), c(1, 0, .5, 0)) # tip 2
  )
  sample_n <- 2
  expect_equal(get_islandMeanFreqM(index_islands, data, sample_n), .5,
               info = "incorrect mean in two tips and two islands, freq .5")
  expect_equal(get_nonislandMeanFreqM(index_nonislands, data, sample_n), .5,
               info = "incorrect mean in two tips and two non-islands, freq .5")
})


test_that("get_SDFreqP structures with equal length islands and non-islands", {
  # One tip and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # 1 tip mean freq p 0.5
  data <- list(c(0, 0, 0), # freq P 0
               c(.5, .5, .5), # freq P 1
               c(0, 0, 0), # freq P 0
               c(.5, .5, .5)) # freq P 1
  exp_SD <- sqrt(4*(0.5)^2/3) 
  sample_n <- 1
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 islands")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 non-islands")
  
  # Two tips and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0), # Tip 1: freq P 0
         c(.5, .5, .5), # Tip 1: freq P 1
         c(0, 0, 0), # Tip 1: freq P 0
         c(.5, .5, .5)), # Tip 1: freq P 1
    list(c(0, 0, 0), # Tip 2: freq P 0
         c(1, .5, 1), # Tip 2: freq P 1/3
         c(.5, 0, .5), # Tip 2: freq P 2/3
         c(.5, .5, .5)) # Tip 2: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tip and 4 islands")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tip and 4 non-islands")
  
  # Two tips and four islands / non-islands (length 1)
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0), # Tip 1: freq P 0
         c(.5), # Tip 1: freq P 1
         c(0), # Tip 1: freq P 0
         c(.5)), # Tip 1: freq P 1
    list(c(0), # Tip 2: freq P 0
         c(1), # Tip 2: freq P 0
         c(.5), # Tip 2: freq P 1
         c(.5)) # Tip 2: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt(4*(0.5)^2/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 islands (length 1)")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 non-islands (length 1)")
  
  
  
  # Two tips with 4 islands and 4 non-islands (8 structures)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5), # Tip 1 island: freq P 1
         c(.5, .5, .5), # Tip 1 non-island: freq P 1
         c(0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5),# Tip 1 island: freq P 1
         c(.5, .5, .5)), # Tip 1 non-island: freq P 1
    list(c(0, 0, 0), # Tip 2 island: freq P 0
         c(0, 0, 0), # Tip 2 non-island: freq P 0
         c(1, .5, 1), # Tip 2 island: freq P 1/3
         c(1, .5, 1), # Tip 2 non-island: freq P 1/3
         c(.5, 0, .5), # Tip 2 island: freq P 2/3
         c(.5, 0, .5), # Tip 2 non-island: freq P 2/3
         c(.5, .5, .5), # Tip 2 island: freq P 1
         c(.5, .5, .5)) # Tip 2 non-island: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures)")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures)")
  
  # Two tips with only 1 island and 7 non-islands 
  index_islands <- c(1)
  index_nonislands <- c(2, 3, 4, 5, 6, 7, 8)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5), # Tip 1 non-island: freq P 1
         c(.5, .5, .5), # Tip 1 non-island: freq P 1
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5),# Tip 1 non-island: freq P 1
         c(.5, .5, .5)), # Tip 1 non-island: freq P 1
    list(c(0, 0, 0), # Tip 2 island: freq P 0
         c(0, 0, 0), # Tip 2 non-island: freq P 0
         c(1, .5, 1), # Tip 2 non-island: freq P 1/3
         c(1, .5, 1), # Tip 2 non-island: freq P 1/3
         c(.5, 0, .5), # Tip 2 non-island: freq P 2/3
         c(.5, 0, .5), # Tip 2 non-island: freq P 2/3
         c(.5, .5, .5), # Tip 2 non-island: freq P 1
         c(.5, .5, .5)) # Tip 2 non-island: freq P 1
  )
  sample_n <- 2
  expect_true(is.na(get_islandSDFreqP(index_islands, data, sample_n)),
              info = "does not return NA when there is only one island structure per tip")
  
})


test_that("get_SDFreqP structures with different length islands and non-islands", {
  # One tip and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # 1 tip mean freq p 0.5
  data <- list(c(0, 0, 0, 0, 0), # freq P 0
               c(.5, .5, .5), # freq P 1
               c(0, 0, 0, 0, 0), # freq P 0
               c(.5, .5, .5)) # freq P 1
  exp_SD <- sqrt(4*(0.5)^2/3) 
  sample_n <- 1
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 islands")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 non-islands")
  
  # Two tips and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1: freq P 0
         c(.5, .5, .5), # Tip 1: freq P 1
         c(0, 0, 0, 0, 0, 0), # Tip 1: freq P 0
         c(.5, .5, .5)), # Tip 1: freq P 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2: freq P 0
         c(1, .5, 1), # Tip 2: freq P 1/3
         c(.5, 0, .5, 1, .5, .5), # Tip 2: freq P 2/3
         c(.5, .5, .5)) # Tip 2: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 islands")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 non-islands")
  
  # Two tips with 4 islands and 4 non-islands (8 structures)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0, 0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5), # Tip 1 island: freq P 1
         c(.5, .5, .5), # Tip 1 non-island: freq P 1
         c(0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5),# Tip 1 island: freq P 1
         c(.5, .5, .5)), # Tip 1 non-island: freq P 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2 island: freq P 0
         c(0, 0, 0, 0, 0, 0), # Tip 2 non-island: freq P 0
         c(1, .5, 1), # Tip 2 island: freq P 1/3
         c(1, .5, 1), # Tip 2 non-island: freq P 1/3
         c(.5, 0, .5), # Tip 2 island: freq P 2/3
         c(.5, 0, .5), # Tip 2 non-island: freq P 2/3
         c(.5, .5, .5), # Tip 2 island: freq P 1
         c(.5, .5, .5)) # Tip 2 non-island: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures)")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures)")
  
  # Two tips with 4 islands and 4 non-islands (8 structures, last one only one site)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq p 0.5, tip 2 mean freq p 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0, 0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5), # Tip 1 island: freq P 1
         c(.5, .5, .5), # Tip 1 non-island: freq P 1
         c(0, 0, 0), # Tip 1 island: freq P 0
         c(0, 0, 0), # Tip 1 non-island: freq P 0
         c(.5, .5, .5),# Tip 1 island: freq P 1
         c(.5)), # Tip 1 non-island: freq P 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2 island: freq P 0
         c(0, 0, 0, 0, 0, 0), # Tip 2 non-island: freq P 0
         c(1, .5, 1), # Tip 2 island: freq P 1/3
         c(1, .5, 1), # Tip 2 non-island: freq P 1/3
         c(.5, 0, .5), # Tip 2 island: freq P 2/3
         c(.5, 0, .5), # Tip 2 non-island: freq P 2/3
         c(.5, .5, .5), # Tip 2 island: freq P 1
         c(.5)) # Tip 2 non-island: freq P 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqP(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures, last one only one site)")
  expect_equal(get_nonislandSDFreqP(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures, last one only one site)")
  
})


test_that("get_SDFreqM structures with equal length islands and non-islands", {
  # One tip and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # 1 tip mean freq m 0.5
  data <- list(c(0, 0, 0), # freq M 0
               c(1, 1, 1), # freq M 1
               c(.5, 0, 0), # freq M 0
               c(1, 1, 1)) # freq M 1
  exp_SD <- sqrt(4*(0.5)^2/3) 
  sample_n <- 1
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 islands")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 non-islands")
  
  # Two tips and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, .5, 0), # Tip 1: freq M 0
         c(1, 1, 1), # Tip 1: freq M 1
         c(.5, 0, 0), # Tip 1: freq M 0
         c(1, 1, 1)), # Tip 1: freq M 1
    list(c(0, 0, 0), # Tip 2: freq M 0
         c(1, .5, 0), # Tip 2: freq M 1/3
         c(1, 1, .5), # Tip 2: freq M 2/3
         c(1, 1, 1)) # Tip 2: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 islands")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 non-islands")
  
  # Two tips and four islands / non-islands (length 1)
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0), # Tip 1: freq M 0
         c(1), # Tip 1: freq M 1
         c(0), # Tip 1: freq M 0
         c(1)), # Tip 1: freq M 1
    list(c(0), # Tip 2: freq M 0
         c(.5), # Tip 2: freq M 0
         c(1), # Tip 2: freq M 1
         c(1)) # Tip 2: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt(4*(0.5)^2/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 islands (length 1)")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 non-islands (length 1)")
  
  
  
  # Two tips with 4 islands and 4 non-islands (8 structures)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, 0, 0), # Tip 1 island: freq M 0
         c(0, .5, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1), # Tip 1 island: freq M 1
         c(1, 1, 1), # Tip 1 non-island: freq M 1
         c(0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1),# Tip 1 island: freq M 1
         c(1, 1, 1)), # Tip 1 non-island: freq M 1
    list(c(0, 0, 0), # Tip 2 island: freq M 0
         c(0, 0, 0), # Tip 2 non-island: freq M 0
         c(1, .5, 0), # Tip 2 island: freq M 1/3
         c(0, .5, 1), # Tip 2 non-island: freq M 1/3
         c(1, 1, .5), # Tip 2 island: freq M 2/3
         c(.5, 1, 1), # Tip 2 non-island: freq M 2/3
         c(1, 1, 1), # Tip 2 island: freq M 1
         c(1, 1, 1)) # Tip 2 non-island: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures)")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures)")
  
  # Two tips with only 1 island and 7 non-islands 
  index_islands <- c(1)
  index_nonislands <- c(2, 3, 4, 5, 6, 7, 8)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1), # Tip 1 non-island: freq M 1
         c(1, 1, 1), # Tip 1 non-island: freq M 1
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1),# Tip 1 non-island: freq M 1
         c(1, 1, 1)), # Tip 1 non-island: freq M 1
    list(c(0, 0, 0), # Tip 2 island: freq M 0
         c(0, 0, 0), # Tip 2 non-island: freq M 0
         c(1, .5, 0), # Tip 2 non-island: freq M 1/3
         c(0, .5, 1), # Tip 2 non-island: freq M 1/3
         c(.5, 1, 1), # Tip 2 non-island: freq M 2/3
         c(1, 1, .5), # Tip 2 non-island: freq M 2/3
         c(1, 1, 1), # Tip 2 non-island: freq M 1
         c(1, 1, 1)) # Tip 2 non-island: freq M 1
  )
  sample_n <- 2
  expect_true(is.na(get_islandSDFreqP(index_islands, data, sample_n)),
              info = "does not return NA when there is only one island structure per tip")
  
})


test_that("get_SDFreqM structures with different length islands and non-islands", {
  # One tip and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # 1 tip mean freq m 0.5
  data <- list(c(0, 0, 0, 0, 0), # freq M 0
               c(1, 1, 1), # freq M 1
               c(0, 0, 0, 0, 0), # freq M 0
               c(1, 1, 1)) # freq M 1
  exp_SD <- sqrt(4*(0.5)^2/3) 
  sample_n <- 1
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 islands")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), exp_SD,
               info = "incorrect SD in one tip and 4 non-islands")
  
  # Two tips and four islands / non-islands
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1: freq M 0
         c(1, 1, 1), # Tip 1: freq M 1
         c(0, 0, 0, 0, 0, 0), # Tip 1: freq M 0
         c(1, 1, 1)), # Tip 1: freq M 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2: freq M 0
         c(1, .5, 0), # Tip 2: freq M 1/3
         c(.5, 0, 1, 1, 1, 1), # Tip 2: freq M 2/3
         c(1, 1, 1)) # Tip 2: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 islands")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect SD in two tips and 4 non-islands")
  
  # Two tips with 4 islands and 4 non-islands (8 structures)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0, 0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1), # Tip 1 island: freq M 1
         c(1, 1, 1), # Tip 1 non-island: freq M 1
         c(0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1),# Tip 1 island: freq M 1
         c(1, 1, 1)), # Tip 1 non-island: freq M 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2 island: freq M 0
         c(0, 0, 0, 0, 0, 0), # Tip 2 non-island: freq M 0
         c(1, .5, 0), # Tip 2 island: freq M 1/3
         c(0, .5, 1), # Tip 2 non-island: freq M 1/3
         c(1, 1, .5), # Tip 2 island: freq M 2/3
         c(.5, 1, 1), # Tip 2 non-island: freq M 2/3
         c(1, 1, 1), # Tip 2 island: freq M 1
         c(1, 1, 1)) # Tip 2 non-island: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures)")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures)")
  
  # Two tips with 4 islands and 4 non-islands (8 structures, last one only one site)
  index_islands <- c(1, 3, 5, 7)
  index_nonislands <- c(2, 4, 6, 8)
  # Tip 1 mean freq m 0.5, tip 2 mean freq m 0.5
  data <- list(
    list(c(0, 0, 0, 0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0, 0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1), # Tip 1 island: freq M 1
         c(1, 1, 1), # Tip 1 non-island: freq M 1
         c(0, 0, 0), # Tip 1 island: freq M 0
         c(0, 0, 0), # Tip 1 non-island: freq M 0
         c(1, 1, 1),# Tip 1 island: freq M 1
         c(1)), # Tip 1 non-island: freq M 1
    list(c(0, 0, 0, 0, 0, 0), # Tip 2 island: freq M 0
         c(0, 0, 0, 0, 0, 0), # Tip 2 non-island: freq M 0
         c(1, .5, 0), # Tip 2 island: freq M 1/3
         c(0, .5, 1), # Tip 2 non-island: freq M 1/3
         c(1, 1, .5), # Tip 2 island: freq M 2/3
         c(.5, 1, 1), # Tip 2 non-island: freq M 2/3
         c(1, 1, 1), # Tip 2 island: freq M 1
         c(1)) # Tip 2 non-island: freq M 1
  )
  exp_SD1 <- sqrt(4*(0.5)^2/3)
  exp_SD2 <- sqrt((2*(0.5)^2+2*((1/3)-0.5)^2)/3)
  mean_SD <- mean(c(exp_SD1, exp_SD2))
  sample_n <- 2
  expect_equal(get_islandSDFreqM(index_islands, data, sample_n), mean_SD,
               info = "incorrect island SD in two tips 4 islands and 4 non-islands (8 structures, last one only one site)")
  expect_equal(get_nonislandSDFreqM(index_nonislands, data, sample_n), mean_SD,
               info = "incorrect non-island SD in two tips 4 islands and 4 non-islands (8 structures, last one only one site)")
  
})



test_that("meanCor", {
  # Expect NA with structures under min_CpG
  
  index_islands <- c(1, 2, 3, 4)
  index_nonislands <- c(1, 2, 3, 4)
  # 1 tip
  data <- list(c(0, 0, 0, 0, .5), 
               c(.5, 1, 1), 
               c(.5, 0, 0, 0, .5), 
               c(1, 1, 1)) 
  expect_true(is.na(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 0, data, sample_n = 1)),
              info = "meanCor_i returns non-NA value when all structures are smaller than 'min_CpG' (single tip)")
  expect_true(is.na(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 0, data, sample_n = 1)),
              info = "meanCor_ni returns non-NA value when all structures are smaller than 'min_CpG' (single tip)")
  # 2 tips
  data <- list(
    list(c(0, 0, 0, 0, 0, .5), 
         c(.5, 1, 1), 
         c(.5, .5, 0, 0, 0, .5), 
         c(.5, 1, 1)), 
    list(c(0, 0, 0, 0, 0, .5), 
         c(.5, 1, 1), 
         c(.5, .5, 0, 0, 0, .5), 
         c(.5, 1, 1))
  )
  expect_true(is.na(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 0, data, sample_n = 2)),
              info = "meanCor_i returns non-NA value when all structures are smaller than 'min_CpG' (2 tips)")
  expect_true(is.na(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 0, data, sample_n = 2)),
              info = "meanCor_ni returns non-NA value when all structures are smaller than 'min_CpG' (2 tips)")
  
  # Expect NA when there is no variation within structures
  # 1 tip
  data <- list(c(0, 0, 0, 0, 0), 
               c(1, 1, 1, 1, 1), 
               c(0, 0, 0, 0, 0), 
               c(1, 1, 1)) 
  expect_true(is.na(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 1)),
              info = "meanCor_i returns non-NA value when there is no variation within structures (single tip)")
  expect_true(is.na(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 1)),
              info = "meanCor_ni returns non-NA value when there is no variation within structures (single tip)")
  
  # 2 tips
  data <- list(
    list(c(0, 0, 0, 0, 0), 
         c(1, 1, 1, 1, 1), 
         c(0, 0, 0, 0, 0), 
         c(1, 1, 1)), 
    list(c(0, 0, 0, 0, 0), 
         c(1, 1, 1, 1, 1), 
         c(0, 0, 0, 0, 0), 
         c(1, 1, 1))
  )
  expect_true(is.na(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 2)),
              info = "meanCor_i returns non-NA value when there is no variation within structures (2 tips)")
  expect_true(is.na(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 2)),
              info = "meanCor_ni returns non-NA value when there is no variation within structures (2 tips)")
  
  
  get_cor <- function(s1, s2){
    cov(s1,s2)/(sd(s1)*sd(s2))
  }
  # Test case shore length 0 equal sizes
  # 1 tip
  data <- list(c(.5, 0, 0, 0, .5), 
               c(.5, 1, 1), 
               c(.5, 0, 0, .5, 1), 
               c(1, 1, 1)) 
  cor1 <- get_cor(data[[1]][1:4], data[[1]][2:5]) ## -0.33
  cor2 <- get_cor(data[[3]][1:4], data[[3]][2:5]) ## 0.301
  expMeanCor <- mean(c(cor1, cor2)) ## -0.0159
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 1), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 0 1 tip)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 1), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 0 1 tip)")
  
  # 2 tips
  data <- list(
    list(c(.5, .5, 0, 0, 0), 
         c(.5, 1, 1), 
         c(.5, .5, 0, .5, .5), 
         c(.5, 1, 1)), 
    list(c(.5, .5, .5, 0, 0), 
         c(.5, 1, 1), 
         c(.5, 0, 0, .5, .5), 
         c(.5, 1, 1))
  )
  cor1_t1 <- get_cor(data[[1]][[1]][1:4], data[[1]][[1]][2:5]) ## 0.5773503
  cor2_t1 <- get_cor(data[[1]][[3]][1:4], data[[1]][[3]][2:5]) ## -0.3333333
  cor1_t2 <- get_cor(data[[2]][[1]][1:4], data[[2]][[1]][2:5]) ## 0.5773503
  cor2_t2 <- get_cor(data[[2]][[3]][1:4], data[[2]][[3]][2:5]) ## 0
  expMeanCor <- mean(c(cor1_t1, cor2_t1, cor1_t2, cor2_t2))
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 2), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 0 2 tip)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 2), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 0 2 tip)")
  
  # Test case shore length 0 non-equal sizes
  # 1 tip
  data <- list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 10 sites
               c(.5, 1, 1, 1, .5), # 5 sites
               c(.5, 0, 0, .5, 1, 1), # 6 sites
               c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, .5, .5)) #12 sites
  cor1 <- get_cor(data[[1]][1:9], data[[1]][2:10]) ## 0.5
  cor2 <- get_cor(data[[2]][1:4], data[[2]][2:5]) ## -0.3333333
  cor3 <- get_cor(data[[3]][1:5], data[[3]][2:6]) ## 0.5976143
  cor4 <- get_cor(data[[4]][1:11], data[[4]][2:12]) ## 0.7370277
  expMeanCor <- mean(c(cor1, cor2, cor3, cor4)) ## 0.3753272
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 1), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 0 1 tip non-equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 1), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 0 1 tip non-equal sizes)")
  
  # 2 tips
  data <- list(
    list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 10 sites
         c(.5, 1, 1, 1, .5), # 5 sites
         c(.5, 0, 0, .5, 1, 1), # 6 sites
         c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, .5, .5)), #12 sites 
    list(c(.5, 0, 0, 0, 0, 0, .5, .5, 1, 1), # 10 sites
         c(.5, .5, 1, .5, .5), # 5 sites
         c(.5, .5, 0, .5, 1, 1), # 6 sites
         c(1, .5, .5, .5, .5, .5, 0, 0, 0, .5, 1, .5)) #12 sites
  )
  cor1_t1 <- get_cor(data[[1]][[1]][1:9], data[[1]][[1]][2:10]) ## 0.5
  cor2_t1 <- get_cor(data[[1]][[2]][1:4], data[[1]][[2]][2:5]) ## -0.3333333
  cor3_t1 <- get_cor(data[[1]][[3]][1:5], data[[1]][[3]][2:6]) ## 0.5976143
  cor4_t1 <- get_cor(data[[1]][[4]][1:11], data[[1]][[4]][2:12]) ## 0.7370277
  cor1_t2 <- get_cor(data[[2]][[1]][1:9], data[[2]][[1]][2:10]) ## 0.7284928
  cor2_t2 <- get_cor(data[[2]][[2]][1:4], data[[2]][[2]][2:5]) ## -0.3333333
  cor3_t2 <- get_cor(data[[2]][[3]][1:5], data[[2]][[3]][2:6]) ## 0.4225771
  cor4_t2 <- get_cor(data[[2]][[4]][1:11], data[[2]][[4]][2:12]) ## 0.4303315
  expMeanCor <- mean(c(cor1_t1, cor2_t1, cor3_t1, cor4_t1, cor1_t2, cor2_t2, cor3_t2, cor4_t2))
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 5, shore_length = 0, data, sample_n = 2), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 0 2 tips non-equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 5, shore_length = 0, data, sample_n = 2), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 0 2 tips non-equal sizes)")
  
  
  # Test case shore length 10 equal sizes
  # 1 tip
  data <- list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 30 sites
               c(.5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 30 sites
               c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1), # 30 sites
               c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5)) # 30 sites
  cor1 <- get_cor(data[[1]][11:19], data[[1]][12:20]) ## 0.5
  cor2 <- get_cor(data[[2]][11:19], data[[2]][12:20]) ## -1.355253e-20
  cor3 <- get_cor(data[[3]][11:19], data[[3]][12:20]) ## 0.2828427
  cor4 <- get_cor(data[[4]][11:19], data[[4]][12:20]) ## 0.7385489
  expMeanCor <- mean(c(cor1, cor2, cor3, cor4)) ## 0.3803479
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 10, data, sample_n = 1), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 10 1 tip equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 10, data, sample_n = 1), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 10 1 tip equal sizes)")
  
  # 2 tips
  data <- list(
    list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 30 sites
         c(.5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 30 sites
         c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1), # 30 sites
         c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5)), # 30 sites
    list(c(.5, 0, 0, .5, .5, .5, 0, 0, .5, 1, .5, 0, 0, 0, 0, .5, .5, 1, 1, 1, .5, 0, 0, 0, .5, .5, 1, 1, 1, 1), # 30 sites
         c(.5, .5, 1, 1, .5, .5, 1, 1, 1, .5, .5, 0, 0, 0, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 30 sites
         c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, .5, .5, 0, 0, .5, 1, 1, 1, 1, .5, 1, .5, 0, 0, .5, 1), # 30 sites
         c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, 1, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, .5, .5, 0, .5)) # 30 sites
  )
  cor1_t1 <- get_cor(data[[1]][[1]][11:19], data[[1]][[1]][12:20]) ## 0.5
  cor2_t1 <- get_cor(data[[1]][[2]][11:19], data[[1]][[2]][12:20]) ## -1.355253e-20
  cor3_t1 <- get_cor(data[[1]][[3]][11:19], data[[1]][[3]][12:20]) ## 0.2828427
  cor4_t1 <- get_cor(data[[1]][[4]][11:19], data[[1]][[4]][12:20]) ## 0.7385489
  cor1_t2 <- get_cor(data[[2]][[1]][11:19], data[[2]][[1]][12:20]) ## 0.7723028
  cor2_t2 <- get_cor(data[[2]][[2]][11:19], data[[2]][[2]][12:20]) ## 0.6666667
  cor3_t2 <- get_cor(data[[2]][[3]][11:19], data[[2]][[3]][12:20]) ## 0.2236068
  cor4_t2 <- get_cor(data[[2]][[4]][11:19], data[[2]][[4]][12:20]) ## 0.7777138
  expMeanCor <- mean(c(cor1_t1, cor2_t1, cor3_t1, cor4_t1, cor1_t2, cor2_t2, cor3_t2, cor4_t2))
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 10, data, sample_n = 2), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 10 2 tips equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 10, data, sample_n = 2), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 10 2 tips equal sizes)")
  
  
  # Test case shore length 10 non-equal sizes
  # 1 tip
  data <- list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 30 sites
               c(.5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 25 sites
               c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1), # 40 sites
               c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, .5, 0, 0, 0, .5)) # 35 sites
  cor1 <- get_cor(data[[1]][11:19], data[[1]][12:20]) ## 0.5
  cor3 <- get_cor(data[[3]][11:29], data[[3]][12:30]) ## 0.2835379
  cor4 <- get_cor(data[[4]][11:24], data[[4]][12:25]) ## 0.6938887
  expMeanCor <- mean(c(cor1, cor3, cor4)) ## 0.4924755
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 10, data, sample_n = 1), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 10 1 tip non-equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 10, data, sample_n = 1), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 10 1 tip non-equal sizes)")
  
  # 2 tips
  data <- list(
    list(c(.5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1, .5, 0, 0, 0, .5, .5, .5, .5, .5, 1), # 30 sites
         c(.5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 25 sites
         c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1), # 40 sites
         c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, 0, 0, 0, .5, .5, 0, 0, 0, .5)), # 35 sites
    list(c(.5, 0, 0, .5, .5, .5, 0, 0, .5, 1, .5, 0, 0, 0, 0, .5, .5, 1, 1, 1, .5, 0, 0, 0, .5, .5, 1, 1, 1, 1), # 30 sites
         c(.5, .5, 1, 1, .5, .5, 1, 1, 1, .5, .5, 0, 0, 0, .5, .5, 1, 1, 1, .5, .5, 1, 1, 1, .5), # 25 sites
         c(.5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, .5, .5, 0, 0, .5, 1, 1, 1, 1, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1, .5, 0, 0, .5, 1), # 40 sites
         c(1, 1, 1, .5, .5, .5, 0, 0, 0, .5, 1, 1, 1, 1, .5, .5, 0, 0, 0, .5, 1, 1, 1, .5, .5, .5, .5, .5, 0, .5, .5, .5, .5, 0, .5)) # 35 sites
  )
  cor1_t1 <- get_cor(data[[1]][[1]][11:19], data[[1]][[1]][12:20]) ## 0.5
  cor3_t1 <- get_cor(data[[1]][[3]][11:29], data[[1]][[3]][12:30]) ## 0.2835379
  cor4_t1 <- get_cor(data[[1]][[4]][11:24], data[[1]][[4]][12:25]) ## 0.6938887
  cor1_t2 <- get_cor(data[[2]][[1]][11:19], data[[2]][[1]][12:20]) ## 0.7723028
  cor3_t2 <- get_cor(data[[2]][[3]][11:29], data[[2]][[3]][12:30]) ## 0.5234868
  cor4_t2 <- get_cor(data[[2]][[4]][11:24], data[[2]][[4]][12:25]) ## 0.7139942
  expMeanCor <- mean(c(cor1_t1, cor3_t1, cor4_t1, cor1_t2, cor3_t2, cor4_t2))
  expect_equal(compute_meanCor_i(index_islands, minN_CpG = 10, shore_length = 10, data, sample_n = 2), expMeanCor,
               info = "MeanCor_i returns non-correct value (shore 10 2 tips non-equal sizes)")
  expect_equal(compute_meanCor_ni(index_nonislands, minN_CpG = 10, shore_length = 10, data, sample_n = 2), expMeanCor,
               info = "MeanCor_ni returns non-correct value (shore 10 2 tips non-equal sizes)")
})

test_that("validate_tree input errors", {
  # No input
  expect_error(validate_tree(),
               info = "function fails to throw an error when no tree is given")
  
  # Tree not character string or phylo class
  tree <- 5
  expect_error(validate_tree(tree),
               info = "function fails to throw error when tree is not character string or phylo object")
  
  # Incorrect newick format
  tree <- "a:1, b:2"
  expect_error(validate_tree(tree),
               info = "function fails to throw error when given tree has incorrect format")
  tree <- "a:1, b:2;"
  expect_error(validate_tree(tree),
               info = "function fails to throw error when given tree has incorrect format")
  # Test error when tree has only one tip
  tree <- "(1:1);"
  expect_error(validate_tree(tree),
               info = "function fails to throw error when given tree has only one tip")
})


test_that("validate_tree returns phylo tree",{
  tree <- "(a:1,b:2);"
  o <- validate_tree(tree)
  expect_equal(class(o), "phylo")
})


test_that("get_cherryDist input control errors", {
  
  # No input
  expect_error(get_cherryDist(),
               info = "function fails to throw an error when no tree is given")
  
})

test_that("get_cherryDist processing of different input types", {
  
  # One cherry, numeric tip labels ordered
  type <- "One cherry, numeric tip labels ordered"
  tree <- "((1:0.1,2:0.1):2,3:3);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, "1",
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, "2",
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, 1,
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, 2,
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, 0.2,
               info = paste("Fails to output correct $dist tip with input type:", type))
  
  # One cherry, numeric tip labels unordered
  type <- "One cherry, numeric tip labels unordered"
  tree <- "((3:0.15,1:0.2):2,2:3);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, "3",
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, "1",
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, 1,
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, 2,
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, 0.35,
               info = paste("Fails to output correct $dist tip with input type:", type))

  
  # One cherry, named tips
  type <- "One cherry, named tips"
  tree <- "((a:0.8,b:0.2):2,c:3);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, "a",
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, "b",
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, 1,
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, 2,
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, 1,
               info = paste("Fails to output correct $dist tip with input type:", type))
  
  
  # One cherry, named tips, cherry to the right of tree
  type <- "One cherry, named tips, cherry to the right of tree"
  tree <- "(a:3,(c:2.1:,b:2.1):1);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, "c",
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, "b",
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, 2,
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, 3,
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, 4.2,
               info = paste("Fails to output correct $dist tip with input type:", type))
  
  
  # Two cherries, numeric tip labels ordered
  type <- "Two cherries, numeric tip labels ordered"
  tree <- "((1:0.1,2:0.1):3,(3:3,4:3):0.2);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, c("1", "3"),
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, c("2", "4"),
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, c(1,3),
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, c(2,4),
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, c(0.2, 6),
               info = paste("Fails to output correct $dist tip with input type:", type))

  
  # Two cherries, numeric tip labels unordered
  type <- "Two cherries, numeric tip labels unordered"
  tree <- "((1:0.15,5:0.2):10,(2:10,10:10):0.35);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, c("1", "2"),
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, c("5", "10"),
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, c(1,3),
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, c(2,4),
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, c(0.35, 20),
               info = paste("Fails to output correct $dist tip with input type:", type))
  
  
  # Two cherries, named tips
  type <- "Two cherries, named tips"
  tree <- "((a:0.15,b:0.2):10,(c:0.1,1d:10):0.35);"
  output <- get_cherryDist(tree)
  expect_equal(output$first_tip_name, c("a", "c"),
               info = paste("Fails to output correct $first_tip_name with input type:", type))
  expect_equal(output$second_tip_name, c("b", "1d"),
               info = paste("Fails to output correct $second_tip_name with input type:", type))
  expect_equal(output$first_tip_index, c(1,3),
               info = paste("Fails to output correct $first_tip_index with input type:", type))
  expect_equal(output$second_tip_index, c(2,4),
               info = paste("Fails to output correct $second_tip_index with input type:", type))
  expect_equal(output$dist, c(0.35, 10.1),
               info = paste("Fails to output correct $dist tip with input type:", type))

})


test_that("validate_data", {
  
  # Different number of structures across tips
  type <- "Different number of structures across tips"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(validate_data(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Tip without structures
  type <- "Tip without structures"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list()) 
  cherryDist <- get_cherryDist(tree)
  expect_error(validate_data(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Data with a structure of length 0
  type <- "Data with a structure of length 0"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), c()), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), c()), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), c())) 
  cherryDist <- get_cherryDist(tree)
  expect_error(validate_data(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Different number of sites in data
  type <- "Different number of sites in data"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(validate_data(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Number of tips in data smaller than given cherry tip indices
  type <- "Number of tips in data smaller than given cherry tip indices"
  tree <- "((1:1,2:1):2,(3:2,5:2):1);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(validate_data(cherryDist, data),
               info = paste(type, "fails to throw error"))
})

test_that("countSites_cherryMethDiff input control errors", {
  
  # Incomplete argument cherryDist
  type <- "Incomplete argument cherryDist"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  cherryDist <- cherryDist[,-1]
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # No cherry info in cherryDist
  type <- "No cherry info in cherryDist"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  cherryDist <- cherryDist[-1,]
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Different number of structures across tips
  type <- "Different number of structures across tips"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Tip without structures
  type <- "Tip without structures"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list()) 
  cherryDist <- get_cherryDist(tree)
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Data with a structure of length 0
  type <- "Data with a structure of length 0"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), c()), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), c()), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), c())) 
  cherryDist <- get_cherryDist(tree)
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Different number of sites in data
  type <- "Different number of sites in data"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
  
  # Number of tips in data smaller than given cherry tip indices
  type <- "Number of tips in data smaller than given cherry tip indices"
  tree <- "((1:1,2:1):2,(3:2,5:2):1);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  expect_error(countSites_cherryMethDiff(cherryDist, data),
               info = paste(type, "fails to throw error"))
})

test_that("countSites_cherryMethDiff correct output with different input types", {
  
  # One cherry numeric labels (left newick)
  type <- "One cherry numeric labels (left newick)"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  cherryDist <- get_cherryDist(tree)
  o <- countSites_cherryMethDiff(cherryDist, data)
  expect_equal(nrow(o), 1, 
               info = paste("Number of output rows does not correspond with input type:", type))
  expect_equal(o$tip_names, "1-2", 
               info = paste("Tip names do not correspond with input type:", type))
  expect_equal(o$tip_indices, "1-2", 
               info = paste("Tip indices do not correspond with input type:", type))
  expect_equal(o$dist, 2, 
               info = paste("Tips distance does not correspond with input type:", type))
  expect_true(all(as.numeric(o[1,-c(1,2,3)]) == c(0, 0, 0, 10, 10, 0)),
              info = paste("Incorrect count number with input type:", type))
  
  # One cherry named labels (right newick)
  type <- "One cherry named labels (right newick)"
  tree <- "(c:5,(b:1.3,a:1.3):2.4);"
  data <- list(
    list(rep(1,10), rep(0,15), rep(1,15)), 
    list(rep(1,10), rep(0.5,15), rep(0,15)), # tip 1 cherry
    list(c(rep(0.5,8), 0, 1), c(rep(0.5,10), 1, 1, 0, 0, 0), rep(0,15))) # tip 2 cherry
  cherryDist <- get_cherryDist(tree)
  o <- countSites_cherryMethDiff(cherryDist, data)
  expect_equal(nrow(o), 1, 
               info = paste("Number of output rows does not correspond with input type:", type))
  expect_equal(o$tip_names, "b-a", 
               info = paste("Tip names do not correspond with input type:", type))
  expect_equal(o$tip_indices, "2-3", 
               info = paste("Tip indices do not correspond with input type:", type))
  expect_equal(o$dist, 2.6, 
               info = paste("Tips distance does not correspond with input type:", type))
  expect_true(all(as.numeric(o[1,-c(1,2,3)]) == c(1, 8, 0, 5, 0, 0)),
              info = paste("Incorrect count number with input type:", type))
  
  # Test case two cherries numeric labels unordered
  tree <- "((25:1.51,24:1.51):2,(38:2.1,45:2.1):1.5);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1, rep(0.5, 5)), c(0.5, 1, rep(0, 8))))
  cherryDist <- get_cherryDist(tree)
  o <- countSites_cherryMethDiff(cherryDist, data)
  expect_equal(nrow(o), 2, 
               info = paste("Number of output rows does not correspond with input type:", type))
  expect_equal(o$tip_names, c("25-24", "38-45"), 
               info = paste("Tip names do not correspond with input type:", type))
  expect_equal(o$tip_indices, c("1-2", "3-4"), 
               info = paste("Tip indices do not correspond with input type:", type))
  expect_equal(o$dist, c(3.02, 4.2), 
               info = paste("Tips distances do not correspond with input type:", type))
  expect_true(all(as.numeric(o[1,-c(1,2,3)]) == c(0, 0, 0, 10, 10, 0)),
              info = "case two cherries incorrect count output for first cherry")
  expect_true(all(as.numeric(o[2,-c(1,2,3)]) == c(5, 5, 0, 5, 1, 1)),
              info = "case two cherries incorrect count output for second cherry")
  
})


test_that("freqSites_cherryMethDiff input control error", {

  # No input
  expect_error(freqSites_cherryMethDiff(),
               info = "function fails to throw an error when no tree is given")
  
  ## Check errors for incorrect tree argument ##
  
  # Tree not character string or phylo class
  tree <- 5
  expect_error(freqSites_cherryMethDiff(tree),
               info = "function fails to throw error when tree is not character string or phylo object")
  
  # Incorrect newick format
  tree <- "a:1, b:2"
  expect_error(freqSites_cherryMethDiff(tree),
               info = "function fails to throw error when given tree has incorrect format")
  tree <- "a:1, b:2;"
  expect_error(freqSites_cherryMethDiff(tree),
               info = "function fails to throw error when given tree has incorrect format")
  # Test error when tree has only one tip
  tree <- "(1:1);"
  expect_error(freqSites_cherryMethDiff(tree),
               info = "function fails to throw error when given tree has only one tip")
  
  
  ## Check errors for incorrect data argument ##
  
  # Different number of structures across tips
  type <- "Different number of structures across tips"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10))) 
  expect_error(freqSites_cherryMethDiff(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Tip without structures
  type <- "Tip without structures"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list()) 
  expect_error(freqSites_cherryMethDiff(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Data with a structure of length 0
  type <- "Data with a structure of length 0"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), c()), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), c()), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), c())) 
  expect_error(freqSites_cherryMethDiff(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Different number of sites in data
  type <- "Different number of sites in data"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  expect_error(freqSites_cherryMethDiff(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Number of tips in data smaller than given cherry tip indices
  type <- "Number of tips in data smaller than given cherry tip indices"
  tree <- "((1:1,2:1):2,(3:2,5:2):1);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  expect_error(freqSites_cherryMethDiff(tree, data),
               info = paste(type, "fails to throw error"))
})


test_that("freqSites_cherryMethDiff output", {
  
  # Test counts are transformed into frequencies correctly
  tree <- "((a:1.5,b:1.5):2,(c:2,d:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,5), rep(1,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
  o <- freqSites_cherryMethDiff(tree = tree, data = data)
  expect_equal(o$tip_names[1], "a-b",
               info = "incorrect tip names in cherry 1")
  expect_equal(o$tip_names[2], "c-d",
               info = "incorrect tip names in cherry 2")
  expect_equal(o$tip_indices[1], "1-2",
               info = "incorrect tip names in cherry 1")
  expect_equal(o$tip_indices[2], "3-4",
               info = "incorrect tip names in cherry 2")
  expect_equal(o$dist[1], 3,
               info = "incorrect tips distance in cherry 1")
  expect_equal(o$dist[2], 4,
               info = "incorrect tips distance in cherry 2")
  column_CpGn <- rep(c(10, 5, 8), each=2)
  f_h_count <- c(0, 0, 0, 5, 8, 0)
  f_h_freq <- f_h_count/column_CpGn
  expect_equal(as.numeric(o[1,4:ncol(o)]), f_h_freq,
               info = "incorrect f and h freqs in cherry 1")
  f_h_count <- c(5, 5, 0, 5, 1, 1)
  f_h_freq <- f_h_count/column_CpGn
  expect_equal(as.numeric(o[2,4:ncol(o)]), f_h_freq,
               info = "incorrect f and h freqs in cherry 2")
})

test_that("get_siteFChange_cherry wrong input", {
  
  # No input
  expect_error(get_siteFChange_cherry(),
               info = "function fails to throw an error when no tree is given")
  
  ## Check errors for incorrect tree argument ##
  
  # Tree not character string or phylo class
  tree <- 5
  expect_error(get_siteFChange_cherry(tree),
               info = "function fails to throw error when tree is not character string or phylo object")
  
  # Incorrect newick format
  tree <- "a:1, b:2"
  expect_error(get_siteFChange_cherry(tree),
               info = "function fails to throw error when given tree has incorrect format")
  tree <- "a:1, b:2;"
  expect_error(get_siteFChange_cherry(tree),
               info = "function fails to throw error when given tree has incorrect format")
  # Test error when tree has only one tip
  tree <- "(1:1);"
  expect_error(get_siteFChange_cherry(tree),
               info = "function fails to throw error when given tree has only one tip")
  
  
  ## Check errors for incorrect data argument ##
  
  # Different number of structures across tips
  type <- "Different number of structures across tips"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10))) 
  expect_error(get_siteFChange_cherry(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Tip without structures
  type <- "Tip without structures"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list()) 
  expect_error(get_siteFChange_cherry(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Data with a structure of length 0
  type <- "Data with a structure of length 0"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), c()), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), c()), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), c())) 
  expect_error(get_siteFChange_cherry(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Different number of sites in data
  type <- "Different number of sites in data"
  tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  expect_error(get_siteFChange_cherry(tree, data),
               info = paste(type, "fails to throw error"))
  
  # Number of tips in data smaller than given cherry tip indices
  type <- "Number of tips in data smaller than given cherry tip indices"
  tree <- "((1:1,2:1):2,(3:2,5:2):1);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,11)), # tip 1 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10)), # tip 2 cherry
    list(rep(1,10), rep(0.5,10), rep(0,10))) 
  expect_error(get_siteFChange_cherry(tree, data),
               info = paste(type, "fails to throw error"))
})

test_that("get_siteFChange_cherry output", {
  
  # Test frequency of sites changing per structure and cherry is computed correctly
  tree <- "((a:1.5,b:1.5):2,(c:2,d:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,5), rep(1,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
  
  o <- get_siteFChange_cherry(tree = tree, data = data)
  expect_equal(o$tip_names[1], "a-b",
               info = "incorrect tip names in cherry 1")
  expect_equal(o$tip_names[2], "c-d",
               info = "incorrect tip names in cherry 2")
  expect_equal(o$tip_indices[1], "1-2",
               info = "incorrect tip names in cherry 1")
  expect_equal(o$tip_indices[2], "3-4",
               info = "incorrect tip names in cherry 2")
  expect_equal(o$dist[1], 3,
               info = "incorrect tips distance in cherry 1")
  expect_equal(o$dist[2], 4,
               info = "incorrect tips distance in cherry 2")
  Fchange_firstCherry <- c(0,1,1)
  expect_equal(as.numeric(o[1,4:ncol(o)]), Fchange_firstCherry,
               info = "incorrect f and h freqs in cherry 1")
  Fchange_secondCherry <- c(1,1,0.25)
  expect_equal(as.numeric(o[2,4:ncol(o)]), Fchange_secondCherry,
               info = "incorrect f and h freqs in cherry 2")
})



