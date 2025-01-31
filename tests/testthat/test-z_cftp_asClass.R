test_that("combiStructureGenerator $get_highest_rate", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  output <- c$get_highest_rate()
  
  # test numerical value
  # test length 1
  # test positive number
  
})

test_that("combiStructureGenerator $get_singleStr_siteNumber and $get_singleStr_number", {
  
  # test 1
  test <- 1
  infoStr <- data.frame(n = c(10, 8, 15),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  expect_equal(c$get_singleStr_number(), length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(c$get_singleStr_siteNumber(), infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
  # test 2
  test <- 2
  infoStr <- data.frame(n = c(25, 78, 1, 16, 88, 40),
                        globalState = c("M", "U", "M", "M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  expect_equal(c$get_singleStr_number(), length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(c$get_singleStr_siteNumber(), infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
})