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

