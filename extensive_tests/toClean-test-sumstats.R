

test_that("get_siteFChange_cherry", {
  tree <- "((1:1.5,2:1.5):2,(3:2,4:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,5), rep(1,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
  o <- get_siteFChange_cherry(tree = tree, data = data, sample_n = 4)
  exp_siteFChange_cherry1 <- c(0, 1, 1)
  exp_siteFChange_cherry2 <- c(1, 1, 2/8)
  expect_equal(as.numeric(o[1,3:ncol(o)]), exp_siteFChange_cherry1,
               info = "incorrect frequency of site FChange cherry 1")
  expect_equal(as.numeric(o[2,3:ncol(o)]), exp_siteFChange_cherry2,
               info = "incorrect frequency of site FChange cherry 2")
})

test_that("SDandMean_SiteFChange_cherry", {
  tree <- "((1:1.5,2:1.5):2,(3:2,4:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,5), rep(1,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
  siteFChange_cherry <- get_siteFChange_cherry(tree = tree, data = data, sample_n = 4)
  
  # Test case all structures are islands
  index_islands <- c(1,2,3)
  index_nonislands <- c()
  o <- SDandMean_SiteFChange_cherry(siteFChange_cherry = siteFChange_cherry, index_islands = index_islands, index_nonislands = index_nonislands)
  exp_mean_I_ch1 <- mean(c(0,1,1))
  exp_mean_I_ch2 <- mean(c(1,1,0.25))
  exp_sd_I_ch1 <- sd(c(0,1,1))
  exp_sd_I_ch2 <- sd(c(1,1,0.25))
  expect_equal(o$mean_I, c(exp_mean_I_ch1, exp_mean_I_ch2),
               info = "returns incorrect mean_I when all str are islands")
  expect_equal(o$sd_I, c(exp_sd_I_ch1, exp_sd_I_ch2),
               info = "returns incorrect sd_I when all str are islands")
  expect_true(all(is.na(o$mean_NI)),
               info = "returns non-NA mean_NI when all str are islands")
  expect_true(all(is.na(o$sd_NI)),
              info = "returns non-NA sd_NI when all str are islands")
  
  # Test case all structures are non-islands
  index_nonislands <- c(1,2,3)
  index_islands <- c()
  o <- SDandMean_SiteFChange_cherry(siteFChange_cherry = siteFChange_cherry, index_islands = index_islands, index_nonislands = index_nonislands)
  exp_mean_NI_ch1 <- mean(c(0,1,1))
  exp_mean_NI_ch2 <- mean(c(1,1,0.25))
  exp_sd_NI_ch1 <- sd(c(0,1,1))
  exp_sd_NI_ch2 <- sd(c(1,1,0.25))
  expect_equal(o$mean_NI, c(exp_mean_NI_ch1, exp_mean_NI_ch2),
               info = "returns incorrect mean_NI when all str are non islands")
  expect_equal(o$sd_NI, c(exp_sd_NI_ch1, exp_sd_NI_ch2),
               info = "returns incorrect sd_NI when all str are non islands")
  expect_true(all(is.na(o$mean_I)),
              info = "returns non-NA mean_I when all str are non islands")
  expect_true(all(is.na(o$sd_I)),
              info = "returns non-NA sd_I when all str are non islands")
  
  # Test case 2 non-islands, 1 island
  index_nonislands <- c(1,3)
  index_islands <- c(2)
  o <- SDandMean_SiteFChange_cherry(siteFChange_cherry = siteFChange_cherry, index_islands = index_islands, index_nonislands = index_nonislands)
  exp_mean_I_ch1 <- mean(c(1))
  exp_mean_NI_ch1 <- mean(c(0,1))
  exp_mean_I_ch2 <- mean(c(1))
  exp_mean_NI_ch2 <- mean(c(1,0.25))
  exp_sd_I_ch1 <- sd(c(1))
  exp_sd_NI_ch1 <- sd(c(0,1))
  exp_sd_I_ch2 <- sd(c(1))
  exp_sd_NI_ch2 <- sd(c(1,0.25))
  expect_equal(o$mean_I, c(exp_mean_I_ch1, exp_mean_I_ch2),
               info = "returns incorrect mean_I with 1 island 2 non-islands")
  expect_equal(o$mean_NI, c(exp_mean_NI_ch1, exp_mean_NI_ch2),
               info = "returns incorrect mean_NI with 1 island 2 non-islands")
  expect_true(all(is.na(o$sd_I)),
              info = "returns non-NA sd_I with 1 island 2 non-islands")
  expect_equal(o$sd_NI, c(exp_sd_NI_ch1, exp_sd_NI_ch2),
               info = "returns incorrect sd_NI with 1 island 2 non-islands")
  
  # Test case 2 islands, 1 non-island
  index_islands <- c(1,3)
  index_nonislands <- c(2)
  o <- SDandMean_SiteFChange_cherry(siteFChange_cherry = siteFChange_cherry, index_islands = index_islands, index_nonislands = index_nonislands)
  exp_mean_NI_ch1 <- mean(c(1))
  exp_mean_I_ch1 <- mean(c(0,1))
  exp_mean_NI_ch2 <- mean(c(1))
  exp_mean_I_ch2 <- mean(c(1,0.25))
  exp_sd_NI_ch1 <- sd(c(1))
  exp_sd_I_ch1 <- sd(c(0,1))
  exp_sd_NI_ch2 <- sd(c(1))
  exp_sd_I_ch2 <- sd(c(1,0.25))
  expect_equal(o$mean_I, c(exp_mean_I_ch1, exp_mean_I_ch2),
               info = "returns incorrect mean_I with 2 islands 1 non-island")
  expect_equal(o$mean_NI, c(exp_mean_NI_ch1, exp_mean_NI_ch2),
               info = "returns incorrect mean_NI with 2 islands 1 non-island")
  expect_equal(o$sd_I, c(exp_sd_I_ch1, exp_sd_I_ch2),
               info = "returns incorrect sd_NI with 2 islands 1 non-island")
  expect_true(all(is.na(o$sd_NI)),
              info = "returns non-NA sd_NI with 2 islands 1 non-island")
  
  
})
