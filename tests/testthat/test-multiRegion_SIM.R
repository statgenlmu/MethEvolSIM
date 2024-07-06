# Function to access private variables and functions
get_private <- function(x) {
  x[['.__enclos_env__']]$private
}

test_that("singleStructureGenerator initialization",{
  # Test my_combiStructure initialization
  single_obj <- singleStructureGenerator$new("U",100)
  expect_true(is.null(get_private(single_obj)$my_combiStructure),
              info = "direct singleStructureGenerator initialization assigns not null $my_combiStructure")

  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState= c("M", "U", "M"))

  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_false(is.null(get_private(combi_obj$get_singleStr(i))$my_combiStructure))
    expect_equal(class(get_private(combi_obj$get_singleStr(i))$my_combiStructure)[1], "combiStructureGenerator",
                 info ="single Str unit initialization by combiStructureGenerator assigns incorrect $my_combiStructure")
    expect_equal(length(get_private(combi_obj)$singleStr[[i]]$get_seq()), 100,
                 info ="single Str unit initialization by combiStructureGenerator generates incorrect sequence length")
    expect_equal(length(get_private(get_private(combi_obj)$singleStr[[1]])$siteR), 100,
                 info ="single Str unit initialization by combiStructureGenerator generates incorrect siteR length")
    expect_equal(length(get_private(get_private(combi_obj)$singleStr[[1]])$eqFreqs), 3,
                 info ="single Str unit initialization by combiStructureGenerator generates incorrect eqFreqs length")
    expect_equal(get_private(combi_obj$get_singleStr(i))$combiStructure_index, i,
                 info ="single Str unit initialization by combiStructureGenerator assigns incorrect combiStructure_index")
  }

})

test_that("combiStructureGenerator initialization", {
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState= c("M", "U", "M"))

  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(length(get_private(combi_obj)$singleStr), nrow(infoStr),
               info = "length of $singleStr not equal to row number in input dataframe")
  expect_equal(length(get_private(combi_obj)$singleStr_globalState), nrow(infoStr),
               info = "length of $singleStr_globalState not equal to row number in input dataframe")
  for (i in 1:nrow(infoStr)){
    expect_equal(class(combi_obj$get_singleStr(i))[1], "singleStructureGenerator",
                 info = "$singleStr elements not of 'singleStructureGenerator' class")
  }
})

test_that("singleStructureGenerator get_seq()", {
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  # Define expected $seq under testing mode
  expected_seq <- c(1, 1, 2, 3, 1, 1, 1, 3, 2, 2, 3, 2, 3)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seq(), expected_seq,
                 info = "returns $seq different from expected under testing mode")
  }

})

test_that("singleStructureGenerator get_seqFirstPos() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seqFirstPos(), 1,
                 info = "function returns $seq first position different from known under testing mode")
  }
  # Under default, non testing mode: unknown $seq
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seqFirstPos(), combi_obj$get_singleStr(i)$get_seq()[1],
                 info = "output not consistent with get_seq()[1]")
  }
})

test_that("singleStructureGenerator get_seqLastPos() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seqLastPos(), 3,
                 info = "function returns $seq last position different from known under testing mode")
  }
  # Under default, non testing mode: unknown $seq
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seqLastPos(), combi_obj$get_singleStr(i)$get_seq()[length(combi_obj$get_singleStr(i)$get_seq())],
                 info = "output not consistent with get_seq()[length(..)] under non testing mode")
  }
})

test_that("singleStructureGenerator get_leftStr_neighbSt() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  expect_error(get_private(combi_obj$get_singleStr(1))$get_leftStr_neighbSt(),
               info = "fails to throw an error when called from first Structure")
  for (i in 2:nrow(infoStr)){
    expect_equal(get_private(combi_obj$get_singleStr(i))$get_leftStr_neighbSt(), 3,
                 info = "function returns left neighbSt different from known under testing mode")
  }
  # Under default, non testing mode: unknown $seq
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 2:nrow(infoStr)){
    expect_equal(get_private(combi_obj$get_singleStr(i))$get_leftStr_neighbSt(), combi_obj$get_singleStr(i-1)$get_seq()[length(combi_obj$get_singleStr(i-1)$get_seq())],
                 info = "output not consistent with get_seq()[length(..)] under non testing mode")
  }
})

test_that("singleStructureGenerator get_rightStr_neighbSt() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  expect_error(get_private(combi_obj$get_singleStr(3))$get_rightStr_neighbSt(),
               info = "fails to throw an error when called from last Structure")
  for (i in 1:(nrow(infoStr)-1)){
    expect_equal(get_private(combi_obj$get_singleStr(i))$get_rightStr_neighbSt(), 1,
                 info = "function returns right neighbSt different from known under testing mode")
  }
  # Under default, non testing mode: unknown $seq
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:(nrow(infoStr)-1)){
    expect_equal(get_private(combi_obj$get_singleStr(i))$get_rightStr_neighbSt(), combi_obj$get_singleStr(i+1)$get_seq()[1],
                 info = "output not consistent with get_seq()[1] under non testing mode")
  }
})

test_that("singleStructureGenerator get_nextStr()", {
  infoStr <- data.frame(n = c(10, 1, 10),
                        globalState= c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(length(get_private(combi_obj$get_singleStr(1))$get_nextStr()$get_seq()), 1)
  expect_equal(length(get_private(combi_obj$get_singleStr(2))$get_nextStr()$get_seq()), 10)
  expect_null(get_private(combi_obj$get_singleStr(3))$get_nextStr())
  single_obj <- singleStructureGenerator$new("U",100)
  expect_null(get_private(single_obj)$get_nextStr())

})


test_that("singleStructureGenerator get_prevStr()", {
  infoStr <- data.frame(n = c(10, 1, 10),
                        globalState= c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(get_private(combi_obj$get_singleStr(1))$get_prevStr())
  expect_equal(length(get_private(combi_obj$get_singleStr(2))$get_prevStr()$get_seq()), 10)
  expect_equal(length(get_private(combi_obj$get_singleStr(3))$get_prevStr()$get_seq()), 1)
  single_obj <- singleStructureGenerator$new("U",100)
  expect_null(get_private(single_obj)$get_prevStr())
})

test_that("singleStructureGenerator init_neighbSt()", {
  # Test cases of singleStructure instances initiated outside combiStructure instance
  obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  expected_neighbSt <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(obj)$neighbSt, expected_neighbSt,
               info = "fails correct neighbSt initialization outside combiStructure instance")

  # Test cases of singleStructure instances initiated from combiStructure instance
  ## Test case 1: initialization of long sequences
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  expected_neighbSt_Str1 <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 4)
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt, expected_neighbSt_Str1,
               info = "assigns incorrect $neighbSt case 1.1")
  expected_neighbSt_Str2 <- c(7, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 4)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt, expected_neighbSt_Str2,
               info = "assigns incorrect $neighbSt case 1.2")
  expected_neighbSt_Str3 <- c(7, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt, expected_neighbSt_Str3,
               info = "assigns incorrect $neighbSt case 1.3")

  # Test case 2: "long" single singleStructure instance initiated from combiStructure instance
  mapNeighbSt_matrix = matrix(c(1L:9L), byrow = TRUE, nrow = 3)
  infoStr <- data.frame(n = c(5),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 2, position 1")
  # position 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 2, position 2")
  # position 3
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 2, position 3")
  # position 4
  leftN <- combi_obj$get_singleStr(1)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 2, position 4")
  # position 5
  leftN <- combi_obj$get_singleStr(1)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 2, position 5")

  # Test case 3: "short" single singleStructure instance initiated from combiStructure instance
  infoStr <- data.frame(n = c(2),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 3, position 1")
  # position 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 3, position 2")

  # Test case 4: "short" singleStructure instances within combiStructure instances with more than 1 singleStructure
  ## Test case 4.1: "short" singleStructure instance is first one
  infoStr <- data.frame(n = c(2, 5, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # singleStr1: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr1: position 1")
  # singleStr1: position 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr1: position 2")
  # singleStr2: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr2: position 1")
  # singleStr2: position 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr2: position 2")
  # singleStr2: position 3
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr2: position 3")
  # singleStr2: position 4
  leftN <- combi_obj$get_singleStr(2)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr2: position 4")
  # singleStr2: position 5
  leftN <- combi_obj$get_singleStr(2)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr2: position 5")
  # singleStr3: position 1
  leftN <- combi_obj$get_singleStr(2)$get_seq()[5]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr3: position 1")
  # singleStr3: position 2
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr3: position 2")
  # singleStr3: position 3
  leftN <- combi_obj$get_singleStr(3)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr3: position 3")
  # singleStr3: position 4
  leftN <- combi_obj$get_singleStr(3)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr3: position 4")
  # singleStr3: position 5
  leftN <- combi_obj$get_singleStr(3)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.1, singleStr3: position 5")

  ## Test case 4.2: "short" singleStructure instance is second one
  infoStr <- data.frame(n = c(5, 2, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # singleStr1: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr1: position 1")
  # singleStr1: position 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr1: position 2")
  # singleStr1: position 3
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr1: position 3")
  # singleStr1: position 4
  leftN <- combi_obj$get_singleStr(1)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr1: position 4")
  # singleStr1: position 5
  leftN <- combi_obj$get_singleStr(1)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr1: position 5")
  # singleStr2: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[5]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr2: position 1")
  # singleStr2: position 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr2: position 2")
  # singleStr3: position 1
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr3: position 1")
  # singleStr3: position 2
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr3: position 2")
  # singleStr3: position 3
  leftN <- combi_obj$get_singleStr(3)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr3: position 3")
  # singleStr3: position 4
  leftN <- combi_obj$get_singleStr(3)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr3: position 4")
  # singleStr3: position 5
  leftN <- combi_obj$get_singleStr(3)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.2, singleStr3: position 5")

  ## Test case 4.3: "short" singleStructure instance is last one
  infoStr <- data.frame(n = c(5, 5, 2),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # singleStr1: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr1: position 1")
  # singleStr1: position 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr1: position 2")
  # singleStr1: position 3
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr1: position 3")
  # singleStr1: position 4
  leftN <- combi_obj$get_singleStr(1)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr1: position 4")
  # singleStr1: position 5
  leftN <- combi_obj$get_singleStr(1)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr1: position 5")
  # singleStr2: position 1
  leftN <- combi_obj$get_singleStr(1)$get_seq()[5]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr2: position 1")
  # singleStr2: position 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr2: position 2")
  # singleStr2: position 3
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr2: position 3")
  # singleStr2: position 4
  leftN <- combi_obj$get_singleStr(2)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[4], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr2: position 4")
  # singleStr2: position 5
  leftN <- combi_obj$get_singleStr(2)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[5], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr2: position 5")
  # singleStr3: position 1
  leftN <- combi_obj$get_singleStr(2)$get_seq()[5]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[1], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr3: position 1")
  # singleStr3: position 2
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "incorrect initialization test case 4.3, singleStr3: position 2")

})

test_that("init_neighbSt() in structure of length 1", {
  mapNeighbSt_matrix = matrix(c(1L:9L), byrow = TRUE, nrow = 3)
  # Test case 1: structure of length 1 in singleStructure instance
  obj <- singleStructureGenerator$new("M", 1)
  expect_equal(get_private(obj)$neighbSt, mapNeighbSt_matrix[obj$get_seq(), obj$get_seq()],
               info ="assigns incorrect $neighbSt in singleStructure instance of length 1")

  # Test case 2: intermediate structure of length 1 in combiStructure instance
  infoStr <- data.frame(n = c(10, 1, 10),
                        globalState= c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # Check lengths according to expected
  expect_equal(length(combi_obj$get_singleStr(1)$get_seq()), 10)
  expect_equal(length(combi_obj$get_singleStr(2)$get_seq()), 1)
  expect_equal(length(combi_obj$get_singleStr(3)$get_seq()), 10)
  # Check neighbSt accordint to expected
  leftN <- combi_obj$get_singleStr(1)$get_seq()[10]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt, exp_neighbSt,
               info = "intermediate singleStr of length 1 gets wrong assignment of neighbSt")

  # Test case 3: first structure of length 1 in combiStructure instance
  infoStr <- data.frame(n = c(1, 10, 10),
                        globalState= c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # Check lengths according to expected
  expect_equal(length(combi_obj$get_singleStr(1)$get_seq()), 1)
  expect_equal(length(combi_obj$get_singleStr(2)$get_seq()), 10)
  expect_equal(length(combi_obj$get_singleStr(3)$get_seq()), 10)
  # Check neighbSt accordint to expected (first position next structure counts as both neighbors)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt, exp_neighbSt,
               info = "first singleStr of length 1 gets wrong assignment of neighbSt")
  # Test case 3.1: first combiStr instance with only singleStr instance of length 1
  infoStr <- data.frame(n = c(1),
                        globalState= c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  leftN <- combi_obj$get_singleStr(1)$get_seq()
  rightN <- combi_obj$get_singleStr(1)$get_seq()
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt, exp_neighbSt,
               info = "3.1. first and only singleStr of length 1 gets wrong assignment of neighbSt")


  # Test case 4: last structure of length 1 in combiStructure instance
  infoStr <- data.frame(n = c(10, 10, 1),
                        globalState= c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # Check lengths according to expected
  expect_equal(length(combi_obj$get_singleStr(1)$get_seq()), 10)
  expect_equal(length(combi_obj$get_singleStr(2)$get_seq()), 10)
  expect_equal(length(combi_obj$get_singleStr(3)$get_seq()), 1)
  # Check neighbSt accordint to expected (last position in previous structure counts as both neighbors)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[length(combi_obj$get_singleStr(2)$get_seq())]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[length(combi_obj$get_singleStr(2)$get_seq())]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt, exp_neighbSt,
               info = "last singleStr of length 1 gets wrong assignment of neighbSt")
})

test_that("singleStructureGenerator get_seq2ndPos() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # expected under seq length not = 1, $seq second position
  for ( i in 1: nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seq2ndPos(), combi_obj$get_singleStr(i)$get_seq()[2],
                 info = "returns $seq second position different from expected for length not = 1")
  }
  # expected under seq length = 1 followed by another structure
  infoStr <- data.frame(n = c(13, 1, 4),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_singleStr(2)$get_seq2ndPos(), combi_obj$get_singleStr(3)$get_seq()[1],
               info ="returns $sec second position different expected under seq length = 1 followed by another structure")
  # expected under seq length = 1 NOT followed by another structure
  infoStr <- data.frame(n = c(13, 5, 1),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_singleStr(3)$get_seq2ndPos())
})


test_that("singleStructureGenerator get_seq2ndButLastPos() from combiStructureGenerator instance", {
  # Under testing mode: known $seq
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # expected under seq length not = 1, $seq second position
  for ( i in 1: nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_seq2ndButLastPos(), combi_obj$get_singleStr(i)$get_seq()[length(combi_obj$get_singleStr(i)$get_seq())-1],
                 info = "returns $seq second but last position different from expected for length not = 1")
  }
  # expected under seq length = 1 preceded by another structure
  infoStr <- data.frame(n = c(13, 1, 4),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_singleStr(2)$get_seq2ndButLastPos(), combi_obj$get_singleStr(1)$get_seq()[length(combi_obj$get_singleStr(1)$get_seq())],
               info ="returns $sec second but lastposition different expected under seq length = 1 preceded by another structure")
  # expected under seq length = 1 NOT preceded by another structure
  infoStr <- data.frame(n = c(1, 14, 4),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_singleStr(1)$get_seq2ndButLastPos())
})

test_that("singleStructureGenerator update_interStr_firstNeighbSt()", {
  single_obj <- singleStructureGenerator$new("U",13, testing = TRUE)
  expect_equal(get_private(single_obj)$neighbSt[1], 1)
  single_obj$update_interStr_firstNeighbSt(3,2)
  expect_equal(get_private(single_obj)$neighbSt[1], 8)
  single_obj$update_interStr_firstNeighbSt(2, NULL)
  expect_equal(get_private(single_obj)$neighbSt[1], 5)
})

test_that("singleStructureGenerator update_interStr_lastNeighbSt()", {
  single_obj <- singleStructureGenerator$new("U",13, testing = TRUE)
  expect_equal(get_private(single_obj)$neighbSt[length(get_private(single_obj)$neighbSt)], 5)
  single_obj$update_interStr_lastNeighbSt(1,3)
  expect_equal(get_private(single_obj)$neighbSt[length(get_private(single_obj)$neighbSt)], 3)
  single_obj$update_interStr_lastNeighbSt(NULL,3)
  expect_equal(get_private(single_obj)$neighbSt[length(get_private(single_obj)$neighbSt)], 9)
})

test_that("singleStructureGenerator update_intraStr_neighbSt()", {
  single_obj <- singleStructureGenerator$new("U",13, testing = TRUE)

  # Test cases: incorrect input
  expect_error(get_private(single_obj)$update_intraStr_neighbSt("i"))
  expect_error(get_private(single_obj)$update_intraStr_neighbSt(c(1,2)))
  expect_error(get_private(single_obj)$update_intraStr_neighbSt(14))
  expect_error(get_private(single_obj)$update_intraStr_neighbSt(0))

  # Test cases: change of methylation state in different positions within singleStructure instances
  if (! "modify_seqPos"%in% names(singleStructureGenerator$public_methods)){
    singleStructureGenerator$set("public", "modify_seqPos", function(position, newState) {
      private$seq[position] <-newState
    })
  }

  ## Structure of length 1
  single_obj <- singleStructureGenerator$new("U", 1)
  single_obj$modify_seqPos(position = 1, newState = 1)
  get_private(single_obj)$update_intraStr_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 1,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 1")
  single_obj$modify_seqPos(position = 1, newState = 2)
  get_private(single_obj)$update_intraStr_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 5,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 2")
  single_obj$modify_seqPos(position = 1, newState = 3)
  get_private(single_obj)$update_intraStr_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 9,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 3")


  ## Position 1
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 1, newState = 2)
  get_private(single_obj)$update_intraStr_neighbSt(1)
  expected_neighbSt <- c(1, 5, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position n
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 13, newState = 1)
  get_private(single_obj)$update_intraStr_neighbSt(13)
  expected_neighbSt <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 7, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position 2
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 2, newState = 3)
  get_private(single_obj)$update_intraStr_neighbSt(2)
  expected_neighbSt <- c(9, 2, 9, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position n-1
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 12, newState = 1)
  get_private(single_obj)$update_intraStr_neighbSt(12)
  expected_neighbSt <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 4, 9, 1)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position between 2 and n-1
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 5, newState = 2)
  get_private(single_obj)$update_intraStr_neighbSt(5)
  expected_neighbSt <- c(1, 2, 3, 5, 7, 4, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)
})

test_that("singleStructureGenerator update_neighbSt()", {
  ##############
  ### Test cases: change of methylation state in different positions within singleStructure instances
  if (! "modify_seqPos"%in% names(singleStructureGenerator$public_methods)){
    singleStructureGenerator$set("public", "modify_seqPos", function(position, newState) {
      private$seq[position] <-newState
    })
  }
  ## Structure of length 1
  single_obj <- singleStructureGenerator$new("U", 1)
  single_obj$modify_seqPos(position = 1, newState = 1)
  get_private(single_obj)$update_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 1,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 1")
  single_obj$modify_seqPos(position = 1, newState = 2)
  get_private(single_obj)$update_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 5,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 2")
  single_obj$modify_seqPos(position = 1, newState = 3)
  get_private(single_obj)$update_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt, 9,
               info = "singleStr of length 1 assigns wrong neighbSt for newState 3")


  ## Position 1
  ### Length to the right > position + 2
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 1, newState = 2)
  get_private(single_obj)$update_neighbSt(1)
  expected_neighbSt <- c(1, 5, 3, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt,
               info = "wrong singleStr update after modification of position 1 with length to the right > position + 2")
  ### length to the right < position + 2 (seq length = 2: only neighbor of position 2 is position 1)
  single_obj <- singleStructureGenerator$new("U", 2)
  single_obj$modify_seqPos(position = 1, newState = 2)
  get_private(single_obj)$update_neighbSt(1)
  expect_equal(get_private(single_obj)$neighbSt[2], 5,
               info = "wrong singleStr update after modification of position 1 with length to the right < position + 2")

  ## Position n
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 13, newState = 1)
  get_private(single_obj)$update_neighbSt(13)
  expected_neighbSt <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 5, 7, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position 2
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 2, newState = 3)
  get_private(single_obj)$update_neighbSt(2)
  expected_neighbSt <- c(9, 2, 9, 4, 7, 1, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position n-1
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 12, newState = 1)
  get_private(single_obj)$update_neighbSt(12)
  expected_neighbSt <- c(1, 2, 3, 4, 7, 1, 3, 2, 8, 6, 4, 9, 1)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ## Position between 2 and n-1
  single_obj <- singleStructureGenerator$new("U", 13, testing = TRUE)
  single_obj$modify_seqPos(position = 5, newState = 2)
  get_private(single_obj)$update_neighbSt(5)
  expected_neighbSt <- c(1, 2, 3, 5, 7, 4, 3, 2, 8, 6, 5, 9, 5)  # Manually calculated based on the testing sequence
  expect_equal(get_private(single_obj)$neighbSt, expected_neighbSt)

  ##########
  #### Test cases: change of methylation state in different positions within combiStructure instances
  ### Case 2: First position
  ## Case 2.1. First position of structure without left neighbouring structure
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 1, newState = 3)
  #debug(get_private(combi_obj$get_singleStr(1))$update_neighbSt)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(1)
  # Update right neighbSt
  ### a) Length to the right > position + 2
  exp_neighbSt <- 8 #modified left neighb 3, right neighb 2
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.1 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(2),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(1)
  exp_neighbSt <- 5 #first position only singleStr counts as both neighbors
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.1 b)")
  ### c) Length to the right < position + 2 with right structure
  mapNeighbSt_matrix = matrix(c(1L:9L), byrow = TRUE, nrow = 3)
  infoStr <- data.frame(n = c(2, 6),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 1, newState = 2)
  leftN <- 2
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(1)
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.1 c)")




  ## Case 2.2. First position of structure with left neighbouring structure with length > 1
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, testing = TRUE)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  # left neighbSt
  exp_neighbSt <- 5 #left neighb 2, modified right neighb 2
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[length(get_private(combi_obj$get_singleStr(1))$neighbSt)], exp_neighbSt,
               info = "assigns wrong left neighbSt case 2.2")
  # right neighbSt
  ### a) Length to the right > position + 2
  exp_neighbSt <- 5 #modified left neighb 2, right neighb 2
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.2 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(13, 2),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.2 b)")
  ### c) Length to the right < position + 2 with right structure
  infoStr <- data.frame(n = c(13, 2, 3),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.2 c)")

  ## Case 2.3. First position of structure with left neighbouring structure with length = 1 and additional left neighbStr
  infoStr <- data.frame(n = c(13, 1, 5),
                       globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt
  leftN <- combi_obj$get_singleStr(1)$get_seq()[length(combi_obj$get_singleStr(1)$get_seq())]
  combi_obj$get_singleStr(3)$modify_seqPos(position = 1, newState = 2)
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(3))$update_neighbSt(1)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 2.3")
  # right neighbSt
  ### a) Length to the right > position + 2
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.3 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(13, 1, 2),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(3)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(3))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.3 b)")
  ### c) Length to the right < position + 2 with right structure
  infoStr <- data.frame(n = c(13, 1, 2, 1),
                        globalState = c("M", "U", "M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(3)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(3))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(3)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(4)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.3 c)")


  ## Case 2.4. First position of structure with left neighbouring structure with length = 1 without additional left neighbStr
  infoStr <- data.frame(n = c(1, 13, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt (there is no left neighbor, so right neighb counts as both neighbors)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  #debug(get_private(combi_obj$get_singleStr(2))$update_neighbSt)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 2.4")
  # right neighbSt
  ### a) Length to the right > position + 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.4 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(1, 2),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.4 b)")
  ### c) Length to the right < position + 2 with right structure
  infoStr <- data.frame(n = c(1, 2, 1),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 1, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(1)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong right neighbSt case 2.4 c)")


  ### Case 3: 2nd position
  ## Case 3.1. Second position of structure without left neighbouring structure
  infoStr <- data.frame(n = c(14, 5),
                        globalState = c("U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt (there is no left neighbor, so right neighb counts as both neighbors)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 3.1")
  # right neighbSt
  ### a) Length to the right > position + 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.1 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(3),
                        globalState = c("U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.1 b)")
  ### c) Length to the right < position + 2 with right structure
  infoStr <- data.frame(n = c(3, 1),
                        globalState = c("U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.1 c)")

  ## Case 3.2. Second position of structure with left neighbouring structure
  infoStr <- data.frame(n = c(1, 13, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 3.2")
  # right neighbSt
  ### a) Length to the right > position + 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.2 a)")
  ### b) Length to the right < position + 2 and no right structure
  infoStr <- data.frame(n = c(1, 3),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.2 b)")
  ### c) Length to the right < position + 2 with right structure
  infoStr <- data.frame(n = c(1, 3, 2),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong right neighbSt case 3.2 c)")


  ### Case 4: 2ndButLast position
  ## Case 4.1. 2ndButLast position of structure without right neighbouring structure
  infoStr <- data.frame(n = c(14, 5),
                        globalState = c("U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 4, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(4)
  # right neighbSt (there is no right neighbor, so left neighb counts as both neighbors)
  leftN <- combi_obj$get_singleStr(2)$get_seq()[4]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[5], exp_neighbSt,
               info = "assigns wrong right neighbSt case 4.1")
  # left neighbSt
  ### a) Length to the left > position - 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[4]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[3], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.1 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(3),
                        globalState = c("U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.1 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 3),
                        globalState = c("U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.1 c)")


  ## Case 4.2. 2ndButLast position of structure with right neighbouring structure
  infoStr <- data.frame(n = c(1, 13, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 12, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(12)
  # right neighbSt
  leftN <- combi_obj$get_singleStr(2)$get_seq()[12]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[13], exp_neighbSt,
               info = "assigns wrong right neighbSt case 4.2")
  # left neighbSt
  ### a) Length to the left > position - 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[10]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[12]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[11], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.2 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(3, 1),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.2 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 3, 1),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 4.2 c)")


  ### Case 5: Last position
  ## Case 5.1. Last position of structure without right neighbouring structure
  infoStr <- data.frame(n = c(1, 13, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt
  ### a) Length to the left > position - 2
  combi_obj$get_singleStr(3)$modify_seqPos(position = 5, newState = 2)
  leftN <- combi_obj$get_singleStr(3)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(3))$update_neighbSt(5)
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[4], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.1 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(2),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.1 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 2),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.1 c)")


  ## Case 5.2. Last position of structure with right neighbouring structure with length > 1
  infoStr <- data.frame(n = c(1, 13, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 13, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(13)
  # right neighbSt
  leftN <- combi_obj$get_singleStr(2)$get_seq()[13]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(3))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong right neighbSt case 5.2")
  # left neighbSt
  ### a) Length to the left > position - 2
  leftN <- combi_obj$get_singleStr(2)$get_seq()[11]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[13]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[12], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.2 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(2, 12),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.2 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 2, 11),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.2 c)")


  ## Case 5.3. Last position of structure with right neighbouring structure with length = 1 and additional right neighbStr
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 13, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(13)
  # right neighbSt
  leftN <- combi_obj$get_singleStr(1)$get_seq()[13]
  rightN <- combi_obj$get_singleStr(3)$get_seq()[1]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong right neighbSt case 5.3")
  # left neighbSt
  ### a) Length to the left > position - 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[11]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[13]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[12], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.3 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(2, 1, 6),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.3 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 2, 1, 5),
                        globalState = c("M", "U", "M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.3 c)")


  ## Case 5.4. Last position of structure with right neighbouring structure with length = 1 without additional right neighbStr
  infoStr <- data.frame(n = c(13, 1),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 13, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(13)
  # right neighbSt
  leftN <- combi_obj$get_singleStr(1)$get_seq()[13]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[13]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong right neighbSt case 5.4")
  # left neighbSt
  ### a) Length to the left > position - 2
  leftN <- combi_obj$get_singleStr(1)$get_seq()[11]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[13]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[12], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.4 a)")
  ### b) Length to the left < position - 2 and no left structure
  infoStr <- data.frame(n = c(2, 1),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(1)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[2]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.4 b)")
  ### c) Length to the left < position - 2 with left structure
  infoStr <- data.frame(n = c(1, 2, 1),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  combi_obj$get_singleStr(2)$modify_seqPos(position = 2, newState = 2)
  get_private(combi_obj$get_singleStr(2))$update_neighbSt(2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(2)$get_seq()[2]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(2))$neighbSt[1], exp_neighbSt,
               info = "assigns wrong left neighbSt case 5.4 c)")


  ### Case 6: Position between 2 and n-1
  infoStr <- data.frame(n = c(13, 1),
                        globalState = c("M", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # left neighbSt
  combi_obj$get_singleStr(1)$modify_seqPos(position = 3, newState = 2)
  leftN <- combi_obj$get_singleStr(1)$get_seq()[1]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[3]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  get_private(combi_obj$get_singleStr(1))$update_neighbSt(3)
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[2], exp_neighbSt,
               info = "assigns wrong left neighbSt case 6")
  # right neighbSt
  leftN <- combi_obj$get_singleStr(1)$get_seq()[3]
  rightN <- combi_obj$get_singleStr(1)$get_seq()[5]
  exp_neighbSt <- mapNeighbSt_matrix[leftN, rightN]
  expect_equal(get_private(combi_obj$get_singleStr(1))$neighbSt[4], exp_neighbSt,
               info = "assigns wrong right neighbSt case 6")

})

test_that("singleStructureGenerator init_Ri_values()",{
  # singleStructure instance
  single_obj <- singleStructureGenerator$new("U", 5)
  iota <- get_private(single_obj)$iota
  expect_equal(mean(get_private(single_obj)$Ri_values), iota,
               info ="mean Ri value does not correspond to iota in isolated singleStructure instance")
  expect_equal(length(get_private(single_obj)$Ri_values), 3,
               info = "incorrect number of Ri values in isolated singleStructure instance")
  # combiStructure instance
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(mean(get_private(combi_obj$get_singleStr(i))$Ri_values), iota,
                 info ="mean Ri value does not correspond to iota in singleStructure instance within combiStructure")
    expect_equal(length(get_private(combi_obj$get_singleStr(i))$Ri_values), 3,
                 info = "incorrect number of Ri values in singleStructure instance within combiStructure")
  }
})

test_that("singleStructureGenerator init_Rc_values()", {
  # singleStructure instance
  single_obj <- singleStructureGenerator$new("U", 5)
  iota <- get_private(single_obj)$iota
  expect_equal(sum(get_private(single_obj)$Rc_values$Rcl, get_private(single_obj)$Rc_values$Rcr), 1-iota,
               info ="sum of Rc values does not correspond to 1 - iota in isolated singleStructure instance")
  expect_equal(length(get_private(single_obj)$Rc_values), 2,
               info ="incorrect number of Rc values in isolated singleStructure instance")
  # combiStructure instance
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(sum(get_private(combi_obj$get_singleStr(i))$Rc_values$Rcl, get_private(combi_obj$get_singleStr(i))$Rc_values$Rcr), 1-iota,
                 info ="mean Ri value does not correspond to iota in singleStructure instance within combiStructure")
    expect_equal(length(get_private(combi_obj$get_singleStr(i))$Rc_values), 2,
                 info = "incorrect number of Ri values in singleStructure instance within combiStructure")
  }
})

test_that("singleStructureGenerator set_Qi()", {
  # singleStructure instance
  ## Check Qi structure
  obj <- singleStructureGenerator$new("U", 10)
  Qi_list <- get_private(obj)$Qi
  expect_true(is.list(Qi_list), info = "does not generate list in isolated singleStructure instance")
  expect_equal(length(Qi_list), 3, info = "does not generate correct number of rate matrices in isolated singleStructure instance")
  expect_true(all(is.matrix(Qi_list[[1]]), is.matrix(Qi_list[[2]]), is.matrix(Qi_list[[3]])), info = "elements are not matrices in isolated singleStructure instance")
  expect_equal(c(length(Qi_list[[1]]), length(Qi_list[[2]]), length(Qi_list[[3]])), c(9, 9, 9), info = "matrices are not of correct length in isolated singleStructure instance")
  ## Check Qi properties (rate matrices)
  validationStates <- listRateMatrix_validation(Qi_list, "set_Qi() test multiRegion_SIM.R in isolated singleStructure instance")
  output <- listMatrices_validationResults(validationStates)
  expect_null(output, info = "validation properties of Qi not met in isolated singleStructure instance")

  # combiStructure instance
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  ## Check Qi structure
  for (i in 1:nrow(infoStr)){
    Qi_list <- get_private(combi_obj$get_singleStr(i))$Qi
    expect_true(is.list(Qi_list), info = "does not generate listin singleStructure instance within combiStructure")
    expect_equal(length(Qi_list), 3, info = "does not generate correct number of rate matrices in singleStructure instance within combiStructure")
    expect_true(all(is.matrix(Qi_list[[1]]), is.matrix(Qi_list[[2]]), is.matrix(Qi_list[[3]])), info = "elements are not matrices in singleStructure instance within combiStructure")
    expect_equal(c(length(Qi_list[[1]]), length(Qi_list[[2]]), length(Qi_list[[3]])), c(9, 9, 9), info = "matrices are not of correct length in singleStructure instance within combiStructure")
  }
  ## Check Qi properties (rate matrices)
  for (i in 1:nrow(infoStr)){
    Qi_list <- get_private(combi_obj$get_singleStr(i))$Qi
    validationStates <- listRateMatrix_validation(Qi_list, "set_Qi() test multiRegion_SIM.R in singleStructure instance within combiStructure")
    output <- listMatrices_validationResults(validationStates)
    expect_null(output, info = "validation properties of Qi not met in singleStructure instance within combiStructure")
  }
})

test_that("singleStructureGenerator set_Qc()", {
  # singleStructure instance
  ## Check Qc structure
  obj <- singleStructureGenerator$new("U", 10)
  Qc_list <- get_private(obj)$Qc
  expect_true(is.list(Qc_list), info = "does not generate list in isolated singleStructure instance")
  expect_equal(length(Qc_list), 9, info = "does not generate correct number of rate matrices in isolated singleStructure instance")
  expect_true(all(is.matrix(Qc_list[[1]]), is.matrix(Qc_list[[2]]), is.matrix(Qc_list[[3]]),
                  is.matrix(Qc_list[[4]]), is.matrix(Qc_list[[5]]), is.matrix(Qc_list[[6]]),
                  is.matrix(Qc_list[[7]]), is.matrix(Qc_list[[8]]), is.matrix(Qc_list[[9]])),
              info = "elements are not matrices in isolated singleStructure instance")
  expect_equal(c(length(Qc_list[[1]]), length(Qc_list[[2]]), length(Qc_list[[3]]),
                 length(Qc_list[[4]]), length(Qc_list[[5]]), length(Qc_list[[6]]),
                 length(Qc_list[[7]]), length(Qc_list[[8]]), length(Qc_list[[9]])),
               rep(9, 9), info = "matrices are not of correct length in isolated singleStructure instance")
  ## Check Qc properties (rate matrices)
  validationStates <- listRateMatrix_validation(Qc_list, "set_Qc() test multiRegion_SIM.R in isolated singleStructure instance")
  output <- listMatrices_validationResults(validationStates)
  expect_null(output, info = "validation properties of Qc not met in isolated singleStructure instance")

  # combiStructure instance
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  ## Check Qc structure
  for (i in 1:nrow(infoStr)){
    Qc_list <- get_private(combi_obj$get_singleStr(i))$Qc
    expect_true(is.list(Qc_list), info = "does not generate list in singleStructure instance within combiStructure")
    expect_equal(length(Qc_list), 9, info = "does not generate correct number of rate matrices in singleStructure instance within combiStructure")
    expect_true(all(is.matrix(Qc_list[[1]]), is.matrix(Qc_list[[2]]), is.matrix(Qc_list[[3]]),
                    is.matrix(Qc_list[[4]]), is.matrix(Qc_list[[5]]), is.matrix(Qc_list[[6]]),
                    is.matrix(Qc_list[[7]]), is.matrix(Qc_list[[8]]), is.matrix(Qc_list[[9]])),
                info = "elements are not matrices in singleStructure instance within combiStructure")
    expect_equal(c(length(Qc_list[[1]]), length(Qc_list[[2]]), length(Qc_list[[3]]),
                   length(Qc_list[[4]]), length(Qc_list[[5]]), length(Qc_list[[6]]),
                   length(Qc_list[[7]]), length(Qc_list[[8]]), length(Qc_list[[9]])),
                 rep(9, 9), info = "matrices are not of correct length in singleStructure instance within combiStructure")
  }
  ## Check Qc properties (rate matrices)
  for (i in 1:nrow(infoStr)){
    Qc_list <- get_private(combi_obj$get_singleStr(i))$Qc
    validationStates <- listRateMatrix_validation(Qc_list, "set_Qc() test multiRegion_SIM.R in singleStructure instance within combiStructure")
    output <- listMatrices_validationResults(validationStates)
    expect_null(output, info = "validation properties of Qc not met in singleStructure instance within combiStructure")
  }
})

test_that("singleStructureGenerator set_Q()", {
  # singleStructure instance
  ## Check Q structure
  obj <- singleStructureGenerator$new("U", 10)
  Q_list <- get_private(obj)$Q
  expect_true(is.list(Q_list), info = "does not generate list in isolated singleStructure instance")
  expect_true(is.list(Q_list[[1]]), info = "does not generate nested list in isolated singleStructure instance")
  expect_true(is.list(Q_list[[2]]), info = "does not generate nested list in isolated singleStructure instance")
  expect_true(is.list(Q_list[[3]]), info = "does not generate nested list in isolated singleStructure instance")
  expect_equal(length(Q_list), 3, info = "does not generate correct number of nested lists in isolated singleStructure instance")
  expect_equal(length(Q_list[[1]]), 9, info = "does not generate correct number of rate matrices in nested list in isolated singleStructure instance")
  expect_equal(length(Q_list[[2]]), 9, info = "does not generate correct number of rate matrices in nested list in isolated singleStructure instance")
  expect_equal(length(Q_list[[3]]), 9, info = "does not generate correct number of rate matrices in nested list in isolated singleStructure instance")
  expect_true(all(is.matrix(Q_list[[1]][[9]]), is.matrix(Q_list[[2]][[5]]), is.matrix(Q_list[[3]][[1]])), info = "elements are not matrices in isolated singleStructure instance")
  ## Check Q properties (rate matrices)
  for (Ri in 1:3){
    validationStates <- listRateMatrix_validation(Q_list[[Ri]], "set_Q() test multiRegion_SIM.R in isolated singleStructure instance")
    output <- listMatrices_validationResults(validationStates)
    expect_null(output, info = "validation properties of Q not met in isolated singleStructure instance")
  }

  # combiStructure instance
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    ## Check Q structure
    Q_list <- get_private(combi_obj$get_singleStr(i))$Q
    expect_true(is.list(Q_list), info = "does not generate list in singleStructure instance within combiStructure")
    expect_true(is.list(Q_list[[1]]), info = "does not generate nested list in singleStructure instance within combiStructure")
    expect_true(is.list(Q_list[[2]]), info = "does not generate nested list in singleStructure instance within combiStructure")
    expect_true(is.list(Q_list[[3]]), info = "does not generate nested list in singleStructure instance within combiStructure")
    expect_equal(length(Q_list), 3, info = "does not generate correct number of nested lists in singleStructure instance within combiStructure")
    expect_equal(length(Q_list[[1]]), 9, info = "does not generate correct number of rate matrices in nested list in singleStructure instance within combiStructure")
    expect_equal(length(Q_list[[2]]), 9, info = "does not generate correct number of rate matrices in nested list in singleStructure instance within combiStructure")
    expect_equal(length(Q_list[[3]]), 9, info = "does not generate correct number of rate matrices in nested list in singleStructure instance within combiStructure")
    expect_true(all(is.matrix(Q_list[[1]][[9]]), is.matrix(Q_list[[2]][[5]]), is.matrix(Q_list[[3]][[1]])), info = "elements are not matrices in singleStructure instance within combiStructure")
    ## Check Q properties (rate matrices)
    for (Ri in 1:3){
      validationStates <- listRateMatrix_validation(Q_list[[Ri]], "set_Q() test multiRegion_SIM.R in singleStructure instance within combiStructure")
      output <- listMatrices_validationResults(validationStates)
      expect_null(output, info = "validation properties of Q not met in singleStructure instance within combiStructure")
    }
  }
})

test_that("singleStructureGenerator initialize_ratetree()", {
  # singleStructure instance
  ### a) with sequence of long length
  # Check that all $ratetree levels sum the same total rate after initialize_ratetree
  obj <- singleStructureGenerator$new("U", 100)
  expect_equal(sd(sapply(get_private(obj)$ratetree, sum)), 0,
               info=paste("Different total rate sums in levels of ratetree after isolated singleStructure initialization a):\n", sapply(get_private(obj)$ratetree, sum)))
  ### b) with sequence of length 1
  obj <- singleStructureGenerator$new("U", 1)
  expect_equal(length(get_private(obj)$ratetree), 1,
               info=paste("Incorrect length of ratetree after isolated singleStructure initialization b):\n", sapply(get_private(obj)$ratetree, sum)))
  expect_true(is.numeric(get_private(obj)$ratetree[[1]]),
              info = "ratetree initialization does not output numeric rate after isolated singleStructure initialization b)")
  ### c) with sequence of length 2
  # Check that all $ratetree levels sum the same total rate after initialize_ratetree
  obj <- singleStructureGenerator$new("U", 2)
  expect_equal(sd(sapply(get_private(obj)$ratetree, sum)), 0,
               info=paste("Different total rate sums in levels of ratetree after isolated singleStructure initialization c):\n", sapply(get_private(obj)$ratetree, sum)))
  ### d) with equilibrium frequencies as all methylated or all unmethylated the total rate should be 0
  obj <- singleStructureGenerator$new(globalState = "U", 100, eqFreqs = c(1,0,0))
  expect_equal(get_private(obj)$ratetree[[1]][1], 0,
               info = paste("Completely unmethylated sequence has total rate different from 0. ratetree[[1]]:", get_private(obj)$ratetree[[1]]))
  obj <- singleStructureGenerator$new(globalState = "M", 100, eqFreqs = c(0,0,1))
  expect_equal(get_private(obj)$ratetree[[1]][1], 0,
               info = paste("Completely methylated sequence has total rate different from 0. ratetree[[1]]:", get_private(obj)$ratetree[[1]]))

  # combiStructure instance
  ### a) with sequences of long length
  # Check that all $ratetree levels sum the same total rate after initialize_ratetree
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(sd(sapply(get_private(combi_obj$get_singleStr(i))$ratetree, sum)), 0,
                 info=paste("Different total rate sums in levels of ratetree after singleStructure instance within combiStructure initialization:\n", sapply(get_private(obj)$ratetree, sum)))
  }
  ### b) with sequence of length 1
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(length(get_private(combi_obj$get_singleStr(2))$ratetree), 1,
               info=paste("Incorrect length of ratetree after singleStructure instance within combiStructure initialization b):\n", sapply(get_private(obj)$ratetree, sum)))
  expect_true(is.numeric(get_private(combi_obj$get_singleStr(2))$ratetree[[1]]),
              info = "ratetree initialization does not output numeric rate after singleStructure instance within combiStructure initialization b)")
  ### c) with sequence of length 2
  # Check that all $ratetree levels sum the same total rate after initialize_ratetree
  infoStr <- data.frame(n = c(13, 2, 4),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(sd(sapply(get_private(combi_obj$get_singleStr(2))$ratetree, sum)), 0,
               info=paste("Different total rate sums in levels of ratetree after singleStructure instance within combiStructure initialization c):\n", sapply(get_private(obj)$ratetree, sum)))
  
  
})

test_that("singleStructureGenerator update_ratetree()", {
  # singleStructure instance
  ### a) with sequences of long length
  obj <- singleStructureGenerator$new("U", 100)
  newrate <- 0.3
  get_private(obj)$update_ratetree(1, newrate)
  expect_equal(get_private(obj)$ratetree[[length(get_private(obj)$ratetree)]][1], newrate,
               info = "fails to update rate value in isolated singleStructure instance a) first change")
  # Check that all $ratetree levels sum the same total rate after update_ratetree
  expect_equal(sd(sapply(get_private(obj)$ratetree, sum)), 0,
               info=paste("Different total rate sums in levels of ratetree after after one change in isolated singleStructure a):\n", sapply(get_private(obj)$ratetree, sum)))
  newrate <- 0.7
  get_private(obj)$update_ratetree(20, 0.7)
  expect_equal(get_private(obj)$ratetree[[length(get_private(obj)$ratetree)]][20], newrate,
               info = "fails to update rate value in isolated singleStructure instance a) second change")
  expect_equal(sd(sapply(get_private(obj)$ratetree, sum)), 0,
               info=paste("Different total rate sums in levels of ratetree after after 2 changes in isolated singleStructure a):\n", sapply(get_private(obj)$ratetree, sum)))
  ### b) with sequence of length 1
  obj <- singleStructureGenerator$new("U", 1)
  newrate <- 0.3
  get_private(obj)$update_ratetree(1, newrate)
  expect_equal(get_private(obj)$ratetree[[length(get_private(obj)$ratetree)]][1], newrate,
               info = "fails to update rate value in isolated singleStructure instance b) first change")
  newrate <- 0.7
  get_private(obj)$update_ratetree(20, 0.7)
  expect_equal(get_private(obj)$ratetree[[length(get_private(obj)$ratetree)]][20], newrate,
               info = "fails to update rate value in isolated singleStructure instance b) second change")


  # combiStructure instance
  ### a) with sequences of long length
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  get_private(combi_obj$get_singleStr(1))$update_ratetree(1, newrate)
  expect_equal(get_private(combi_obj$get_singleStr(1))$ratetree[[length(get_private(combi_obj$get_singleStr(1))$ratetree)]][1],
               newrate, info = "fails to update rate value in singleStructure instance within combiStructure")
  for (i in 1:nrow(infoStr)){
    expect_equal(sd(sapply(get_private(combi_obj$get_singleStr(i))$ratetree, sum)), 0,
                 info=paste("Different total rate sums in levels of ratetree after singleStructure instance within combiStructure update:\n", sapply(get_private(obj)$ratetree, sum)))
  }
  ### b) with sequence of length 1
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  newrate <- 0.3
  get_private(combi_obj$get_singleStr(2))$update_ratetree(1, newrate)
  expect_equal(get_private(combi_obj$get_singleStr(2))$ratetree[[1]], newrate,
               info = "fails to update rate value in singleStructure instance within combiStructure b) first change")
  newrate <- 0.7
  get_private(combi_obj$get_singleStr(2))$update_ratetree(1, 0.7)
  expect_equal(get_private(combi_obj$get_singleStr(2))$ratetree[[1]], newrate,
               info = "fails to update rate value in singleStructure instance within combiStructure b) second change")
})

test_that("singleStructureGenerator()",{
  obj <- singleStructureGenerator$new("U", 100)
  # Check chosen seqposition is correct
  # with long$seq length
  v <- get_private(obj)$choose_random_seqpos(testing=TRUE)
  index <- v[1]
  sampled_r <- v[2]
  rateTree <- get_private(obj)$ratetree
  rates <- rateTree[[length(rateTree)]]
  expect_true(index == 1 | sum(rates[1:(index-1)]) < sampled_r,
              info = paste("Chosen index too large in long isolated singleStr instance.\n Sum of previous rates:", sum(rates[1:(index-1)]), "Sampled index:", sampled_r))
  expect_true(sampled_r <= sum(rates[1:index]),
              info = paste("Chosen index too small in long isolated singleStr instance.\n Sum of rates from 1 to index:", sum(rates[index:length(rates)]), "Sampled index:", sampled_r))
  # with $seq length of 1
  obj <- singleStructureGenerator$new("U", 1)
  v <- get_private(obj)$choose_random_seqpos(testing=TRUE)
  index <- v[1]
  sampled_r <- v[2]
  rateTree <- get_private(obj)$ratetree
  rate <- rateTree[[length(rateTree)]]
  expect_equal(index, 1,
              info = paste("Chosen index is not only $seq index in isolated singleStr instance of length 1.\n Sum of rates from 1 to index:", sum(rates[index:length(rates)]), "Sampled index:", sampled_r))
})

test_that("singleStructureGenerator choose_number_of_changes()",{
  obj <- singleStructureGenerator$new("U", 100)
  # Check chosen number of changes is integer
  c <- get_private(obj)$choose_number_of_changes(dt = 0.01)
  expect_true(is.integer(c))
  ## If sequence is totally unmethylated or methylated the total rate of change is 0, no changes should be sampled
  obj <- singleStructureGenerator$new(globalState = "U", n = 100, eqFreqs = c(1,0,0))
  dt_possibilities <- seq(from=0.01, to  = 10, by =0.1)
  for(dt in 1:length(dt_possibilities)){
    expect_equal(get_private(obj)$choose_number_of_changes(dt = dt), 0,
                 info = paste("sampled number of changes non-0 for sequence with rate 0"))
  }
  
})

test_that("combiStructureGenerator set_singleStr() and copy()",{
  if (! "modify_seqPos"%in% names(singleStructureGenerator$public_methods)){
    singleStructureGenerator$set("public", "modify_seqPos", function(position, newState) {
      private$seq[position] <-newState
    })
  }
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combiObj <- combiStructureGenerator$new(infoStr)
  cloned_combiObj <- combiObj$copy()
  original_data <- get_private(combiObj$get_singleStr(1))$seq
  copied_data <- get_private(cloned_combiObj$get_singleStr(1))$seq
  expect_equal(original_data, copied_data, info="fails to copy data")

  ## Modify value in cloned_combiObj and check it does not modify combiObj
  cloned_combiObj$get_singleStr(1)$modify_seqPos(1,2)
  changed_value <- get_private(cloned_combiObj$get_singleStr(1))$seq[1]
  expect_equal(changed_value, 2, info = "changed singleStructure instance does not get the given value")
  unchanged_value <- get_private(combiObj$get_singleStr(1))$seq[1]
  expect_equal(unchanged_value, original_data[1], info = "unchanged singleStructure within combiStructure instance changes")

})

test_that("singleStructureGenerator SSE_evol()", {
  ################# singleStructure instance ##
  obj <- singleStructureGenerator$new("U", 100)
  original_obj <- obj$clone()
  # Testing mode: Returns SSE_evol info with dataframe
  number_SSEs <- 100
  dtime = 0.01
  test_info <- list()
  for (i in 1:number_SSEs){
    test_info[[i]] <- obj$SSE_evol(dt = dtime, testing = TRUE)
  }
  # Test 1: Dataframe contains as many rows as events sampled
  for (i in 1:number_SSEs){
    expect_equal(test_info[[i]]$event_number, nrow(test_info[[i]]$SSE_evolInfo),
                 info = paste("Test 1 isolated singleStr failed for SSE number:", i))
  }
  # Check that in positions with changes..
  if(any(sapply(test_info, function(x) x$event_number != 0))){
    indices <- which(sapply(test_info, function(x) x$event_number != 0))
    filtered_dataframes <- do.call(rbind, lapply(indices, function(i) {
      test_info[[i]]$SSE_evolInfo
    }))
    # Test 2: Old sequence equals old state info and new sequence equals new state info
    for (i in filtered_dataframes[,"position"]){
      expect_equal(filtered_dataframes[filtered_dataframes$position==i,"old_St"][1], get_private(original_obj)$seq[i],
                   info = "Test 2 isolated singleStr instance: oldSt of first event in $seq position does not correspond to original sequence")
      expect_equal(filtered_dataframes[filtered_dataframes$position==i,"new_St"][length(filtered_dataframes[filtered_dataframes$position==i,"new_St"])],get_private(obj)$seq[i],
                   info = "Test 2 isolated singleStr instance: newSt of last event in $seq position does not correspond to evolved sequence")
    }
  }
  # Test that by default (testing = FALSE) the function returns nothing
  result <- obj$SSE_evol(dt = 0.01)
  expect_null(result, info="output not null in isolated singleStr instance")

  #################### combiStructure instance ##
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  original_combi_obj <- combi_obj$copy()

  # Testing mode: Returns SSE_evol info with dataframe
  number_SSEs <- 100
  dtime = 0.01
  test_info <- list()
  for (i in 1:number_SSEs){
    test_info[[i]] <- combi_obj$get_singleStr(1)$SSE_evol(dt = dtime, testing = TRUE)
  }
  # Test 1: Dataframe contains as many rows as events sampled
  for (i in 1:number_SSEs){
    expect_equal(test_info[[i]]$event_number, nrow(test_info[[i]]$SSE_evolInfo),
                 info = paste("Test 1 singleStr within combiStr failed for SSE number:", i))
  }
  # Check that in positions with changes..
  if(any(sapply(test_info, function(x) x$event_number != 0))){
    indices <- which(sapply(test_info, function(x) x$event_number != 0))
    filtered_dataframes <- do.call(rbind, lapply(indices, function(i) {
      test_info[[i]]$SSE_evolInfo
    }))
    # Test 2: Old sequence equals old state info and new sequence equals new state info
    for (i in filtered_dataframes[,"position"]){
      expect_equal(filtered_dataframes[filtered_dataframes$position==i,"old_St"][1], get_private(original_combi_obj$get_singleStr(1))$seq[i],
                   info = "Test 2 isolated singleStr instance:oldSt of first event in $seq position does not correspond to original sequence")
      expect_equal(filtered_dataframes[filtered_dataframes$position==i,"new_St"][length(filtered_dataframes[filtered_dataframes$position==i,"new_St"])],get_private(combi_obj$get_singleStr(1))$seq[i],
                   info = "Test 2 isolated singleStr instance: newSt of last event in $seq position does not correspond to evolved sequence")
    }
  }
  # Test that by default (testing = FALSE) the function returns nothing
  result <- combi_obj$get_singleStr(1)$SSE_evol(dt = 0.01)
  expect_null(result, info = "output not null in singleStr instance within combiStr")


})

test_that("combiStructureGenerator SSE_evol()",{
  # Test 1: combiStructure instance with number of singleStr > 1
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  test_info <- get_private(combi_obj)$SSE_evol(dt = 0.01, testing=TRUE)
  expect_equal(length(test_info$evol_order), combi_obj$get_singleStr_number(), info = "Test1: length of singleStr evol_order different from singleStr number")
  expect_equal(length(test_info$testing_info), combi_obj$get_singleStr_number(), info = "Test1: length of SSE testing info different from singleStr number")

  # Test 1: combiStructure instance with number of singleStr > 1
  infoStr <- data.frame(n = c(13),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  test_info <- get_private(combi_obj)$SSE_evol(dt = 0.01, testing=TRUE)
  expect_equal(length(test_info$evol_order), combi_obj$get_singleStr_number(), info = "Test2: length of singleStr evol_order different from singleStr number")
  expect_equal(length(test_info$testing_info), combi_obj$get_singleStr_number(), info = "Test2: length of SSE testing info different from singleStr number")

  # Check it does not produce output with testing = FALSE
  test_info <- get_private(combi_obj)$SSE_evol(dt = 0.01, testing=FALSE)
  expect_null(test_info, info = "outputs info with testing = FALSE")
})

test_that("combiStructureGenerator interval_evol()",{
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  # Testing mode:
  #Interval_length is bigger than dt with remainder
  output <- get_private(combi_obj)$interval_evol(interval_length = 1.005, dt = 0.01, testing=T)
  expect_equal(output[1], paste("Number of SSE:", 100, ". Interval length:", 0.01))
  expect_equal(strsplit(output[2], " ")[[1]][4], "remainder")

  # Interval_length is smaller than dt
  output <- get_private(combi_obj)$interval_evol(interval_length = 0.005, dt = 0.01, testing=T)
  expect_equal(output, paste("SSE evolving in shorter interval_length than dt. Interval length:", 0.005))

  # Interval_length is bigger than dt and no remainder
  output <- get_private(combi_obj)$interval_evol(interval_length = 1, dt = 0.01, testing=T)
  expect_equal(output[1], paste("Number of SSE: 100 . Interval length:", 0.01))
  expect_equal(output[2], "No remainder")

  # Test that by default (testing = FALSE) the function returns nothing
  output <- get_private(combi_obj)$interval_evol(interval_length = 1.005, dt = 0.01)
  expect_null(output)
})

test_that("singleStructureGenerator IWE_evol()", {
  ######## singleStructure instance ##
  ### a) with long $seq length
  obj <- singleStructureGenerator$new("U", 100)
  original_obj <- obj$clone()

  # Testing mode: Returns list with IWE info
  number_IWEs <- 100
  test_info <- list()
  for (i in 1:number_IWEs){
    test_info[[i]] <- obj$IWE_evol(testing = TRUE)
  }
  # Each IWE outputs a list
  expect_true(all(sapply(test_info, is.list)), info ="Not all IWEs with testing=TRUE output list in IWE_evol of isolated singleStructure instance a)")
  # Function to compare Q matrices returns TRUE if not identical and FALSE if identical
  compare_matrices <- function(mat1, mat2) {
    return(!identical(mat1, mat2))
  }
  for(i in 1:number_IWEs){
    expect_true(is.matrix(test_info[[i]]$Mk), info = "testing output $Mk does not return matrix in IWE_evol of isolated singleStructure instance a)")
    if (test_info[[i]]$eqFreqsChange){
      expect_equal(length(test_info[[i]]), 11,
                   info = "IWE output length does not correspond with $eqFreqsChange TRUE in IWE_evol of isolated singleStructure instance a)")
      # Test case matching
      if (test_info[[i]]$IWE_case == "Case 1. u bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[1]>test_info[[i]]$old_eqFreqs[1], info = "fails test case matching u bigger in IWE_evol of isolated singleStructure instance a)")
      }
      if (test_info[[i]]$IWE_case == "Case 1. p bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[2]>test_info[[i]]$old_eqFreqs[2], info = "fails test case matching p bigger in IWE_evol of isolated singleStructure instance a)")
      }
      if (test_info[[i]]$IWE_case == "Case 1. m bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[3]>test_info[[i]]$old_eqFreqs[3], info = "fails test case matching m bigger in IWE_evol of isolated singleStructure instance a)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. u smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[1]<test_info[[i]]$old_eqFreqs[1], info = "fails test case matching u smaller in IWE_evol of isolated singleStructure instance a)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. p smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[2]<test_info[[i]]$old_eqFreqs[2], info = "fails test case matching p smaller in IWE_evol of isolated singleStructure instance a)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. m smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[3]<test_info[[i]]$old_eqFreqs[3], info = "fails test case matching m smaller in IWE_evol of isolated singleStructure instance a)")
      }
      if (i >=2){
        for (Ri in 1:3){
          for (neighbCase in 1:9){
            expect_true(compare_matrices(test_info[[i]]$new_Q[[Ri]][[neighbCase]], test_info[[i-1]]$new_Q[[Ri]][[neighbCase]]),
                        info = "$new_Q not updated with $eqFreqsChange TRUE in IWE_evol of isolated singleStructure instance a)")
          }
        }
      }

    } else {
      expect_false(test_info[[i]]$eqFreqsChange)
      expect_equal(length(test_info[[i]]), 3,
                   info = "IWE output length does not correspond with $eqFreqsChange FALSE in IWE_evol of isolated singleStructure instance a)")
    }
    if (i>=2){
      expect_equal(test_info[[i]]$old_eqFreqs, test_info[[i-1]]$new_eqFreqs,
                   info = "old_eqFreqs after IWE does not correspond with new_eqFreqs in previous IWE in IWE_evol of isolated singleStructure instance a)")
      expect_equal(test_info[[i]]$old_Q, test_info[[i-1]]$new_Q,
                   info = "old_Q after IWE does not correspond with new_Q in previous IWE in IWE_evol of isolated singleStructure instance a)")
    }

  }
  expect_equal(test_info[[1]]$old_eqFreqs, get_private(original_obj)$eqFreqs,
               info = "old_eqFreqs after first IWE does not correspond with eqFreqs original object in IWE_evol of isolated singleStructure instance a)")
  expect_equal(test_info[[number_IWEs]]$new_eqFreqs, get_private(obj)$eqFreqs,
               info = "old_eqFreqs after last IWE does not correspond with eqFreqs modified object in IWE_evol of isolated singleStructure instance a)")

  expect_equal(test_info[[1]]$old_Q, get_private(original_obj)$Q,
               info = "old_Q after first IWE does not correspond with $Q original object in IWE_evol of isolated singleStructure instance a)")
  expect_equal(test_info[[number_IWEs]]$new_Q, get_private(obj)$Q,
               info = "new_Q after last IWE does not correspond with $Q modified object in IWE_evol of isolated singleStructure instance a)")

  obj <- singleStructureGenerator$new("U", 100)
  original_obj <- obj$clone()

  test_info <- obj$IWE_evol(testing = TRUE)

  if (test_info$eqFreqsChange){
    if(!is.null(test_info$changedPos)){
      for (i in test_info$changedPos){
        expect_false(get_private(obj)$seq[i]==get_private(original_obj)$seq[i],
                     info = "$seq[changedPos] not changed with $eqFreqsChange T in IWE_evol of isolated singleStructure instance a)")
        if (i >= 2){
          expect_false(get_private(obj)$neighbSt[i-1]==get_private(original_obj)$neighbSt[i-1],
                       info = "neighbSt not updated in IWE_evol of isolated singleStructure instance a)")
        }
        if (i<= 99){
          expect_false(get_private(obj)$neighbSt[i+1]==get_private(original_obj)$neighbSt[i+1],
                       info = "neighbSt not updated in IWE_evol of isolated singleStructure instance a)")
        }
        for(j in max(i-1, 1):min(i+1, length(get_private(obj)$seq))){
          rate_obj <- abs(get_private(obj)$Q[[get_private(obj)$siteR[j]]][[get_private(obj)$neighbSt[j]]][get_private(obj)$seq[j], get_private(obj)$seq[j]])
          rate_oriObj <- abs(get_private(original_obj)$Q[[get_private(original_obj)$siteR[j]]][[get_private(original_obj)$neighbSt[j]]][get_private(original_obj)$seq[j], get_private(original_obj)$seq[j]])
          if (rate_obj !=rate_oriObj){
            expect_false(get_private(obj)$ratetree[[length(get_private(obj)$ratetree)]][j]==get_private(original_obj)$ratetree[[length(get_private(original_obj)$ratetree)]][j],
                         info = "$ratetree not updated for changedPos and/or neighbours with different change rate in IWE_evol of isolated singleStructure instance a)")
          }
        } #While distances are neglected, there can be a change that leads to a position with neighbSt 2 to neighbSt 4
        # and that does not translate in a different rate because those Qc are equal
      }
    }
  }

  # Test that by default (testing = FALSE) the function returns nothing
  result <- obj$IWE_evol()
  expect_null(result)

  ######## singleStructure instance ##
  ### b) with $seq of length = 1
  obj <- singleStructureGenerator$new("U", 1)
  original_obj <- obj$clone()

  # Testing mode: Returns list with IWE info
  number_IWEs <- 100
  test_info <- list()
  for (i in 1:number_IWEs){
    test_info[[i]] <- obj$IWE_evol(testing = TRUE)
  }

  # Each IWE outputs a list
  expect_true(all(sapply(test_info, is.list)), info ="Not all IWEs with testing=TRUE output list in IWE_evol of isolated singleStructure instance b)")
  for(i in 1:number_IWEs){
    expect_true(is.matrix(test_info[[i]]$Mk), info = "testing output $Mk does not return matrix in IWE_evol of isolated singleStructure instance b)")
    if (test_info[[i]]$eqFreqsChange){
      expect_equal(length(test_info[[i]]), 11,
                   info = "IWE output length does not correspond with $eqFreqsChange TRUE in IWE_evol of isolated singleStructure instance b)")
      # Test case matching
      if (test_info[[i]]$IWE_case == "Case 1. u bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[1]>test_info[[i]]$old_eqFreqs[1], info = "fails test case matching u bigger in IWE_evol of isolated singleStructure instance b)")
      }
      if (test_info[[i]]$IWE_case == "Case 1. p bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[2]>test_info[[i]]$old_eqFreqs[2], info = "fails test case matching p bigger in IWE_evol of isolated singleStructure instance b)")
      }
      if (test_info[[i]]$IWE_case == "Case 1. m bigger"){
        expect_true(test_info[[i]]$new_eqFreqs[3]>test_info[[i]]$old_eqFreqs[3], info = "fails test case matching m bigger in IWE_evol of isolated singleStructure instance b)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. u smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[1]<test_info[[i]]$old_eqFreqs[1], info = "fails test case matching u smaller in IWE_evol of isolated singleStructure instance b)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. p smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[2]<test_info[[i]]$old_eqFreqs[2], info = "fails test case matching p smaller in IWE_evol of isolated singleStructure instance b)")
      }
      if (test_info[[i]]$IWE_case == "Case 2. m smaller"){
        expect_true(test_info[[i]]$new_eqFreqs[3]<test_info[[i]]$old_eqFreqs[3], info = "fails test case matching m smaller in IWE_evol of isolated singleStructure instance b)")
      }
      if (i >=2){
        for (Ri in 1:3){
          for (neighbCase in 1:9){
            expect_true(compare_matrices(test_info[[i]]$new_Q[[Ri]][[neighbCase]], test_info[[i-1]]$new_Q[[Ri]][[neighbCase]]),
                        info = "$new_Q not updated with $eqFreqsChange TRUE in IWE_evol of isolated singleStructure instance b)")
          }
        }
      }

    } else {
      expect_false(test_info[[i]]$eqFreqsChange)
      expect_equal(length(test_info[[i]]), 3,
                   info = "IWE output length does not correspond with $eqFreqsChange FALSE in IWE_evol of isolated singleStructure instance b)")
    }
    if (i>=2){
      expect_equal(test_info[[i]]$old_eqFreqs, test_info[[i-1]]$new_eqFreqs,
                   info = "old_eqFreqs after IWE does not correspond with new_eqFreqs in previous IWE in IWE_evol of isolated singleStructure instance b)")
      expect_equal(test_info[[i]]$old_Q, test_info[[i-1]]$new_Q,
                   info = "old_Q after IWE does not correspond with new_Q in previous IWE in IWE_evol of isolated singleStructure instance b)")
    }

  }
  expect_equal(test_info[[1]]$old_eqFreqs, get_private(original_obj)$eqFreqs,
               info = "old_eqFreqs after first IWE does not correspond with eqFreqs original object in IWE_evol of isolated singleStructure instance b)")
  expect_equal(test_info[[number_IWEs]]$new_eqFreqs, get_private(obj)$eqFreqs,
               info = "old_eqFreqs after last IWE does not correspond with eqFreqs modified object in IWE_evol of isolated singleStructure instance b)")

  expect_equal(test_info[[1]]$old_Q, get_private(original_obj)$Q,
               info = "old_Q after first IWE does not correspond with $Q original object in IWE_evol of isolated singleStructure instance b)")
  expect_equal(test_info[[number_IWEs]]$new_Q, get_private(obj)$Q,
               info = "new_Q after last IWE does not correspond with $Q modified object in IWE_evol of isolated singleStructure instance b)")


})


test_that("combiStructureGenerator get_island_number()",{
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_island_number(), 1)
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("U", "U", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_island_number(), 3)
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_island_number(), 0)
})

test_that("combiStructureGenerator get_island_index()",{
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_island_index(), 2)
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("U", "U", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_equal(combi_obj$get_island_index(), c(1,2,3))
})

test_that("combiStructureGenerator $set_IWE_rate()",{
  # Test 1
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  mu <- get_private(combi_obj)$mu
  island_number <- 1
  expected <- mu * island_number
  expect_equal(get_private(combi_obj)$IWE_rate, expected, info ="test 1 fails")
  # Test 2
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("U", "U", "U"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  island_number <- 3
  expected <- mu * island_number
  expect_equal(get_private(combi_obj)$IWE_rate, expected, info ="test 2 fails")
  # Test 3
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  island_number <- 0
  expected <- mu * island_number
  expect_equal(get_private(combi_obj)$IWE_rate, expected, info ="test 3 fails")

})

test_that("combiStructureGenerator $IWE_events",{
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_IWE_events(), info = "$IWE_events not initialized as NULL")
  combi_obj$set_IWE_events("booh")
  expect_equal(combi_obj$get_IWE_events(), "booh", info = "$get_IWE_events() does not return info assigned with $set_IWE_events()")
})

test_that("combiStructureGenerator $branch_evol()", {
  # Testing mode:
  # No islands: Simulation without IWEs
  infoStr <- data.frame(n = c(100, 100, 100, 100, 100, 100, 100),
                        globalState= c("M", "M", "M", "M", "M", "M", "M"))

  obj <- combiStructureGenerator$new(infoStr)
  branchLength <- 1
  output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
  expect_false(output$IWE_event, info ="without islands it samples IWE_events")
  expect_equal(obj$get_IWE_events(), "Simulation without IWE events.")

  # Long branch length to favor IWE events
  infoStr <- data.frame(n = c(100, 100, 100, 100, 100, 100, 100),
                        globalState= c("M", "U", "M", "U", "U", "U", "U"))

  obj <- combiStructureGenerator$new(infoStr)
  branchLength <- 10
  output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
  if(output$IWE_event){
    expect_true(all(output$IWE_times < branchLength),
                info = "not all IWE_times < branch_length when IWE_event TRUE")
    expect_true(length(output$SSE_intervals) == length(output$IWE_times)+1,
                info = "number of intervals not equal to number of IWE_t + 1 when IWE_event TRUE")
    expect_equal(sum(output$SSE_intervals), branchLength,
                 info = "sum of SSE_interval lengths not equal to branch length when IWE_event TRUE")
    expect_equal(length(output$infoIWE), length(output$IWE_times),
                 info = "infoIWE length not equal to IWE_times length when IWE_event TRUE")
    expect_equal(obj$get_IWE_events()$times, output$IWE_times,
                 info = "$get_IWE_events()$times not equal to testing output IWE_times when IWE_event TRUE")
    expect_equal(obj$get_IWE_events()$islands, output$islands,
                 info = "$get_IWE_events()$islands not equal to testing output IWE_times when IWE_event TRUE")
  } else {
    expect_equal(length(output$IWE_times), 1,
                 info = "length of IWE times not 1 when IWE_event FALSE")
    expect_true(output$IWE_times>branchLength,
                info = "IWE_t is not bigger than branch_length when IWE_event FALSE")
    expect_false(obj$get_IWE_events(),
                 info = "$IWE_events not FALSE when IWE_event FALSE")
  }

  # Short branch length to favor no IWE events
  infoStr <- data.frame(n = c(100, 100, 100, 100, 100, 100, 100),
                        globalState= c("M", "U", "M", "U", "U", "U", "U"))

  obj <- combiStructureGenerator$new(infoStr)
  branchLength <- 1
  output <- obj$branch_evol(branch_length = branchLength, dt = 0.01, testing=T)
  if(output$IWE_event){
    expect_true(all(output$IWE_times < branchLength),
                info = "not all IWE_times < branch_length when IWE_event TRUE")
    expect_true(length(output$SSE_intervals) == length(output$IWE_times)+1,
                info = "number of intervals not equal to number of IWE_t + 1 when IWE_event TRUE")
    expect_equal(sum(output$SSE_intervals), branchLength,
                 info = "sum of SSE_interval lengths not equal to branch length when IWE_event TRUE")
    expect_equal(length(output$infoIWE), length(output$IWE_times),
                 info = "infoIWE length not equal to IWE_times length when IWE_event TRUE")
    expect_equal(obj$get_IWE_events()$times, output$IWE_times,
                 info = "$get_IWE_events()$times not equal to testing output IWE_times when IWE_event TRUE")
    expect_equal(obj$get_IWE_events()$islands, output$islands,
                 info = "$get_IWE_events()$islands not equal to testing output IWE_times when IWE_event TRUE")
  } else {
    expect_equal(length(output$IWE_times), 1,
                 info = "length of IWE times not 1 when IWE_event FALSE")
    expect_true(output$IWE_times>branchLength,
                info = "IWE_t is not bigger than branch_length when IWE_event FALSE")
    expect_false(obj$get_IWE_events(),
                 info = "$IWE_events not FALSE when IWE_event FALSE")
  }

  # Test that by default (testing = FALSE) the function returns nothing
  output <- obj$branch_evol(branch_length = branchLength, dt = 0.01)
  expect_null(output)
})

test_that("combiStructureGenerator $set_name() and $get_name()", {
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_name(), info = "$name is not itialized as null")
  combi_obj$set_name("name")
  expect_equal(combi_obj$get_name(), "name", info = "$set_name() argument does not correspond with $get_name() output")
})

test_that("combiStructureGenerator $set_own_index() and $get_own_index()", {
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_own_index(), info = "$own_index is not itialized as null")
  combi_obj$set_own_index(1)
  expect_equal(combi_obj$get_own_index(), 1, info = "$set_own_index() argument does not correspond with $get_own_index() output")
})

test_that("combiStructureGenerator $set_offspring_index() and $get_offspring_index()", {
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_offspring_index(), info = "$offspring_index is not itialized as null")
  combi_obj$set_offspring_index(1)
  expect_equal(combi_obj$get_offspring_index(), 1, info = "$set_offspring_index() argument does not correspond with $get_offspring_index() output")
})

test_that("combiStructureGenerator $add_offspring_index()", {
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_offspring_index(), info = "$offspring_index is not itialized as null")
  combi_obj$add_offspring_index(1)
  expect_equal(combi_obj$get_offspring_index(), 1, info = "first offspring index")
  combi_obj$add_offspring_index(2)
  expect_equal(combi_obj$get_offspring_index(), c(1,2), info = "second offspring index")
})

test_that("combiStructureGenerator $set_parent_index() and $get_parent_index()", {
  infoStr <- data.frame(n = c(13, 1, 5),
                        globalState = c("M", "M", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$get_parent_index(), info = "$parent_index is not itialized as null")
  combi_obj$set_parent_index(1)
  expect_equal(combi_obj$get_parent_index(), 1, info = "$set_parent_index() argument does not correspond with $get_parent_index() output")
})

test_that("split_newick correctly splits Newick tree", {

  # Test case 1: Test with a simple tree
  tree <- "(a:1, c:2, (d:3.7, e:4):5);"
  expected_output <- data.frame(
    "unit" = c("a", "c", "(d:3.7, e:4)"),
    "brlen" = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 1: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 1: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 1: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 1: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 1: Output should match the expected result.")

  # Test case 2: Test with a subtree
  tree <- "(d:3.7, e:4)"
  expected_output <- data.frame(
    "unit" = c("d", "e"),
    "brlen" = c(3.7, 4),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 2: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 2: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 2: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 2: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 2: Output should match the expected result.")

  # Test case 3: Test with a more complex tree
  tree <- "((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
  expected_output <- data.frame(
    "unit" = c("(a:1, b:1)", "c", "(d:3.7, (e:4, f:1):3)"),
    "brlen" = c(2, 2, 5),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 3: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 3: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 3: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 3: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 3: Output should match the expected result.")

  # Test case 4: Test with complex tree tip names
  tree <- "(weird_name:1, c:2, (d:3.7, 25name_with_numbers0:4):5);"
  expected_output <- data.frame(
    "unit" = c("weird_name", "c", "(d:3.7, 25name_with_numbers0:4)"),
    "brlen" = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 4: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 4: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 4: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 4: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 4: Output should match the expected result.")

  # Test case 5: Test with non integer branch lengths
  tree <- "(a:0.169, c:2, (d:3.782, e:4):5.01);"
  expected_output <- data.frame(
    "unit" = c("a", "c", "(d:3.782, e:4)"),
    "brlen" = c(0.169, 2, 5.01),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 5: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 5: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 5: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 5: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 5: Output should match the expected result.")

  # Test case 6: Test with a tree with different spacing pattern
  tree <- "( a:1, c: 2, (d:3.7, e :4):5);"
  expected_output <- data.frame(
    "unit" = c("a", "c", "(d:3.7, e :4)"),
    "brlen" = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  expect_output <- split_newick(tree)
  expect_true(is.data.frame(expect_output),
              info = "Test case 6: Output should be a data frame.")
  expect_equal(ncol(expect_output), 2,
               info = "Test case 6: Output should have 2 columns.")
  expect_true(is.character(expect_output$unit),
              info = "Test case 6: 'unit' column should be character type.")
  expect_true(is.numeric(expect_output$brlen),
              info = "Test case 6: 'brlen' column should be numeric type.")
  expect_equal(expect_output, expected_output,
               info = "Test case 6: Output should match the expected result.")
})

test_that("treeMultiRegionSimulator", {
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState = c("M", "U", "M"))
  message <- capture.output(treeData <- treeMultiRegionSimulator$new(infoStr, tree = "(a:1, c:2, (d:3.7, e:4):5);"), type = "message")
  expect_equal(message, "Simulating data at root and letting it evolve along given tree:  (a:1, c:2, (d:3.7, e:4):5);")
  expect_equal(length(treeData$Branch), 6, info = "Generates incorrect number of branches")
  expect_equal(length(treeData$Branch), 6, info = "Generates incorrect number of branches")
  expect_null(treeData$Branch[[1]]$get_parent_index(), info = "Data at root has not null parent_index")
  expect_equal(treeData$Branch[[2]]$get_parent_index(), 1, info ="Assigns incorrect parent_index")
  expect_equal(treeData$Branch[[3]]$get_parent_index(), 1, info ="Assigns incorrect parent_index")
  expect_equal(treeData$Branch[[4]]$get_parent_index(), 1, info ="Assigns incorrect parent_index")
  expect_equal(treeData$Branch[[1]]$get_offspring_index(), c(2,3,4), info = "Assigns incorrect offspring_index")
  expect_null(treeData$Branch[[2]]$get_offspring_index(), info = "Assigns incorrect offspring_index")
  expect_null(treeData$Branch[[3]]$get_offspring_index(), info = "Assigns incorrect offspring_index")
  expect_equal(treeData$Branch[[4]]$get_offspring_index(), c(5,6), info = "Assigns incorrect offspring_index")
  expect_null(treeData$Branch[[5]]$get_offspring_index(), info = "Assigns incorrect offspring_index")
  expect_null(treeData$Branch[[6]]$get_offspring_index(), info = "Assigns incorrect offspring_index")
  expect_null(treeData$Branch[[1]]$get_name(), info = "Assigns incorrect name")
  expect_equal(treeData$Branch[[2]]$get_name(), "a", info = "Assigns incorrect name")
  expect_equal(treeData$Branch[[3]]$get_name(), "c", info = "Assigns incorrect name")
  expect_null(treeData$Branch[[4]]$get_name(), info = "Assigns incorrect name")
  expect_equal(treeData$Branch[[5]]$get_name(), "d", info = "Assigns incorrect name")
  expect_equal(treeData$Branch[[6]]$get_name(), "e", info = "Assigns incorrect name")
  expect_equal(length(treeData$branchLength), 6, info = "Generates incorrect number of branchLength")
  expect_equal(treeData$branchLength[1], NA_real_, info = "Assigns incorrect branchLength")
  expect_equal(treeData$branchLength[2], 1, info = "Assigns incorrect branchLength")
  expect_equal(treeData$branchLength[3], 2, info = "Assigns incorrect branchLength")
  expect_equal(treeData$branchLength[4], 5, info = "Assigns incorrect branchLength")
  expect_equal(treeData$branchLength[5], 3.7, info = "Assigns incorrect branchLength")
  expect_equal(treeData$branchLength[6], 4, info = "Assigns incorrect branchLength")
  expect_null(treeData$Branch[[1]]$get_IWE_events(), info = "Data at root has not null IWE_events")
  for (branch in 2:6){
    expect_false(is.null(treeData$Branch[[branch]]$get_IWE_events()), info = "Branches should have $IWE_events FALSE/vector of times")
  }
  for (branch in 1:6){
    expect_equal(treeData$Branch[[branch]]$get_singleStr_number(), 3, info = "Branches have incorrect number of singleStructures")
    for (str in 1:3){
      expect_equal(length(treeData$Branch[[branch]]$get_singleStr(str)$get_seq()), 100, info ="Incorrect sequence length")
    }
  }
})

test_that("customized params", {
  # From singleStructure initialization
  params <- get_parameterValues()
  params$alpha_pI <- 0.25
  obj <- singleStructureGenerator$new("U", 10, params = params)
  expect_equal(obj$get_alpha_pI(), 0.25, info = "fails to initiate singleStructure with new parameter value")

  # From combiStructure initialization
  params <- get_parameterValues()
  params$alpha_Ri <- 0.3
  params$mu <- 0.05
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr, params)
  expect_equal(combi_obj$get_singleStr(1)$get_alpha_Ri(), 0.3,
               info = "fails to initiate first structure initialized within combiStructure with new parameter value")
  expect_equal(combi_obj$get_singleStr(2)$get_alpha_Ri(), 0.3,
               info = "fails to initiate second structure initialized within combiStructure with new parameter value")
  expect_equal(combi_obj$get_singleStr(3)$get_alpha_Ri(), 0.3,
               info = "fails to initiate third structure initialized within combiStructure with new parameter value")
  expect_equal(combi_obj$get_mu(), 0.05,
               info ="fails to initiate combiStructure with new parameter value")

  # From treeMultiRegionSimulator initialization
  params <- get_parameterValues()
  params$mu <- 0.002
  params$iota <- 0.1
  message <- capture.output(treeData <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1, c:2, (d:3.7, e:4):5);", params = params), type = "message")
  expect_equal(message, "Simulating data at root and letting it evolve along given tree:  (a:1, c:2, (d:3.7, e:4):5);")
  for (br in 1:6){
    expect_equal(treeData$Branch[[br]]$get_mu(), 0.002,
                 info = paste("fails to initiate new parameter value in combiStructure initialized from treeMultiRegionSimulator, branch number: ", i))
    for(str in 1:3){
      expect_equal(treeData$Branch[[br]]$get_singleStr(str)$get_iota(), 0.1,
                   info = paste("fails to initiate new parameter value in singleStr", str, "of branch", br, "initialized from treeMultiRegionSimulator"))
    }
  }

})

test_that("fixed eqFreqs",{
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0.1, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
  combi_obj <- combiStructureGenerator$new(infoStr = infoStr)
  for (i in 1:nrow(infoStr)){
    expect_equal(combi_obj$get_singleStr(i)$get_eqFreqs(), c(infoStr$u_eqFreq[i], infoStr$p_eqFreq[i], infoStr$m_eqFreq[i]),
                 info = "does not assign correct eqFreqs when given. Generated from combiStructureGenerator")
  }
  silence <- capture.output(tree_obj <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1, c:2, (d:3.7, e:4):5);"), type = "message")
  for (i in 1:nrow(infoStr)){
    expect_equal(tree_obj$Branch[[1]]$get_singleStr(i)$get_eqFreqs(), c(infoStr$u_eqFreq[i], infoStr$p_eqFreq[i], infoStr$m_eqFreq[i]),
                 info = "does not assign correct eqFreqs when given. Generated from combiStructureGenerator")
  }

  # Incorrect frequencies or missing values
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(NA, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
  expect_error(capture.output(combiStructureGenerator$new(infoStr = infoStr), type = "message"),
               info = "fails throwing an error when eqFreqs have missing values. Generated from combiStructureGenerator")
  expect_error(capture.output(treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1, c:2, (d:3.7, e:4):5);"), type = "message"),
               info = "fails throwing an error when eqFreqs have missing values. Generated from treeMultiRegionSimulator")
  infoStr <- data.frame(n = c(100, 100, 100),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
  expect_error(combiStructureGenerator$new(infoStr = infoStr),
               info = "fails throwing an error when eqFreqs wrong frequencies. Generated from combiStructureGenerator")
  expect_error(capture.output(treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1, c:2, (d:3.7, e:4):5);"), type = "message"),
               info = "fails throwing an error when eqFreqs wrong frequencies. Generated from treeMultiRegionSimulator")
})




