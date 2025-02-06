

test_that("get_cherryDist", {
  # Set tree
  newick_tree <- "((1:1,2:1):1,3:2);"
  
  # Test wrong input
  expect_error(get_cherryDist(),
               info = "function fails to throw an error when no tree/sample_n is given")
  expect_error(get_cherryDist(sample_n = 2),
               info = "function fails to throw an error when no tree is given")
  expect_error(get_cherryDist(tree = newick_tree),
               info = "function fails to throw an error when no sample_n is given")
  
  
  # Test error when tree has only one tip
  newick_tree <- "(1:1);"
  expect_error(get_cherryDist(tree = newick_tree, sample_n = 1),
               info = "fails to throw an error when number of tips is < 2")
  
  # Test conversion from newick tree to ape tree
  newick_tree <- "((1:1,2:1):1,3:2);"
  expect_true(class(get_cherryDist(newick_tree, sample_n = 3, testing = T)$tree) == "phylo",
              info = "fails to convert newick tree into ape's class phylo")
  
  # Test case single cherry
  expect_true(all(get_cherryDist(newick_tree, sample_n = 3)[[1]] == c(1,2,2)),
              info = "incorrect output for newick tree with single cherry")
  ape_tree <- ape::read.tree(text = newick_tree)
  expect_true(all(get_cherryDist(ape_tree, sample_n = 3)[[1]] == c(1,2,2)),
              info = "incorrect output for ape tree with single cherry")
  
  # Test case 2 cherries
  newick_tree <- "((1:1.5,2:1.5):2,(3:2,4:2):1.5);"
  ape_tree <- ape::read.tree(text = newick_tree)
  expect_true(all(get_cherryDist(newick_tree, sample_n = 4)[[1]] == c(1,2,3)),
              info = "incorrect output for newick tree with two cherries (1st cherry)")
  expect_true(all(get_cherryDist(newick_tree, sample_n = 4)[[2]] == c(3,4,4)),
              info = "incorrect output for newick tree with two cherries (2nd cherry)")
  expect_true(all(get_cherryDist(ape_tree, sample_n = 4)[[1]] == c(1,2,3)),
              info = "incorrect output for ape tree with two cherries (1st cherry)")
  expect_true(all(get_cherryDist(ape_tree, sample_n = 4)[[2]] == c(3,4,4)),
              info = "incorrect output for ape tree with two cherries (2nd cherry)")
  
})


test_that("countSites_cherryMethDiff", {
  # Test case one cherry
  newick_tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)))
  cherries <- get_cherryDist(newick_tree, sample_n = 3)
  o <- countSites_cherryMethDiff(cherries, data)
  expect_equal(nrow(o), 1, info = "case one cherry does not have one row")
  expect_true(all(as.numeric(o[1,-c(1,2)]) == c(0, 0, 0, 10, 10, 0)),
              info = "case one cherry incorrect count output")
  
  # Test case two cherries
  newick_tree <- "((1:1.5,2:1.5):2,(3:2,4:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1, rep(0.5, 5)), c(0.5, 1, rep(0, 8))))
  cherries <- get_cherryDist(newick_tree, sample_n = 4)
  o <- countSites_cherryMethDiff(cherries, data)
  expect_equal(nrow(o), 2, info = "case two cherries does not have one row")
  expect_true(all(as.numeric(o[1,-c(1,2)]) == c(0, 0, 0, 10, 10, 0)),
              info = "case two cherries incorrect count output for first cherry")
  expect_true(all(as.numeric(o[2,-c(1,2)]) == c(5, 5, 0, 5, 1, 1)),
              info = "case two cherries incorrect count output for second cherry")
  
})

test_that("get_siteFChange_cherry_perMethDiffType", {
  
  # Test error when tree has only one tip
  newick_tree <- "(1:1);"
  data <- list(rep(1,10), rep(0,10), rep(1,10))
  expect_error(get_siteFChange_cherry_perMethDiffType(tree = newick_tree, data = data, sample_n = 1),
               info = "fails to throw an error when number of tips is < 2")
  
  # Test error when arguments are missing
  newick_tree <- "((1:1,2:1):1,3:2);"
  data <- list(
    list(rep(1,10), rep(0,10), rep(1,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)),
    list(rep(1,10), rep(0.5,10), rep(0,10)))
  expect_error(get_siteFChange_cherry_perMethDiffType(data = data, sample_n = 3),
               info = "fails to throw an error when argument 'tree' is not given")
  expect_error(get_siteFChange_cherry_perMethDiffType(tree = newick_tree, sample_n = 3),
               info = "fails to throw an error when argument 'data' is not given")
  expect_error(get_siteFChange_cherry_perMethDiffType(tree = newick_tree, data = data),
               info = "fails to throw an error when argument 'sample_n' is not given")
  
  # Test counts are transformed into frequencies correctly
  tree <- "((1:1.5,2:1.5):2,(3:2,4:2):1.5);"
  data <- list(
    list(rep(1,10), rep(0,5), rep(1,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(rep(1,10), rep(0.5,5), rep(0,8)),
    list(c(rep(0,5), rep(0.5, 5)), c(0, 0, 1, 1, 1), c(0.5, 1, rep(0, 6))))
  o <- get_siteFChange_cherry_perMethDiffType(tree = tree, data = data, sample_n = 4)
  expect_equal(o$tips[1], "1-2",
               info = "incorrect tips id in cherry 1")
  expect_equal(o$tips[2], "3-4",
               info = "incorrect tips id in cherry 2")
  expect_equal(o$dist[1], 3,
               info = "incorrect tips distance in cherry 1")
  expect_equal(o$dist[2], 4,
               info = "incorrect tips distance in cherry 2")
  column_CpGn <- rep(c(10, 5, 8), each=2)
  f_h_count <- c(0, 0, 0, 5, 8, 0)
  f_h_freq <- f_h_count/column_CpGn
  expect_equal(as.numeric(o[1,3:ncol(o)]), f_h_freq,
               info = "incorrect f and h freqs in cherry 1")
  f_h_count <- c(5, 5, 0, 5, 1, 1)
  f_h_freq <- f_h_count/column_CpGn
  expect_equal(as.numeric(o[2,3:ncol(o)]), f_h_freq,
               info = "incorrect f and h freqs in cherry 2")
})

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
