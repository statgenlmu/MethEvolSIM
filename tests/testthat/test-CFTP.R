test_that("singleStructureGenerator $get_Qi()",{
  
  # Initialize singleStructureGenerator instance
  single_obj <- singleStructureGenerator$new("U",10)
  
  # Get list of 3 SSEi rate matrices
  SSEi_3mat <- single_obj$get_Qi()
  
  expect_equal(length(SSEi_3mat), 3,
               info = "method with null arguments does not return 3 objects")
  expect_true(is.matrix(SSEi_3mat[[1]]),
              info = "method with null arguments does not return a list of matrices")
  expect_true(is.matrix(SSEi_3mat[[2]]),
              info = "method with null arguments does not return a list of matrices")
  expect_true(is.matrix(SSEi_3mat[[3]]),
              info = "method with null arguments does not return a list of matrices")
  
  # Get rate for serveral cases
  SSEi_111 <- single_obj$get_Qi(siteR = 1, oldSt = 1, newSt = 1)
  SSEi_213 <- single_obj$get_Qi(siteR = 2, oldSt = 1, newSt = 3)
  SSEi_323 <- single_obj$get_Qi(siteR = 3, oldSt = 2, newSt = 3)
  
  expect_equal(SSEi_111, SSEi_3mat[[1]][1,1],
               info = "method with arguments does not return rate as in matrices")
  expect_equal(SSEi_213, SSEi_3mat[[2]][1,3],
               info = "method with arguments does not return rate as in matrices")
  expect_equal(SSEi_323, SSEi_3mat[[3]][2,3],
               info = "method with arguments does not return rate as in matrices")
})

test_that("singleStructureGenerator $get_leftneighbSt() and $get_rightneighbSt()",{
  
  # Initialize singleStructureGenerator instance
  single_obj <- singleStructureGenerator$new("U",10)
  
  # Get left neighbSt for site 3 and right neighbSt for site 1: that means state of site 2
  expect_equal(single_obj$get_leftneighbSt(3), single_obj$get_neighbSt(2),
               info = "neighbSt obtained with get_leftneighbSt is not correct")
  expect_equal(single_obj$get_rightneighbSt(1), single_obj$get_neighbSt(2),
               info = "neighbSt obtained with get_rightneighbSt is not correct")
  
})

test_that("singleStructureGenerator $cftp_all_equal()",{
  
  # Initialize singleStructureGenerator instance
  single_obj <- singleStructureGenerator$new("U",10)
  
  # Expect error when function is called without giving a value for the 'state' argument
  expect_error(single_obj$cftp_all_equal(),
               info = "method fails to throw error when no value for 'state' argument is given")
  
  # Expect error when value of argument 'state' is not correct
  expect_error(single_obj$cftp_all_equal(state=1),
               info = "method fails to throw error when 'state' value is not correct")
  
  # Expect NULL output when state is correct but testing is (as default) FALSE
  expect_null(single_obj$cftp_all_equal(state="U"),
              info = "whith testing = FALSE method generates output")
  
  # Check generated sequence for both possible states
  expect_equal(single_obj$cftp_all_equal(state = "U", testing = TRUE), rep(1, 10),
               info = "method fails to assign sequence of 'u' states")
  expect_equal(single_obj$cftp_all_equal(state = "M", testing = TRUE), rep(3, 10),
               info = "method fails to assign sequence of 'm' states")
  
})

test_that("singleStructureGenerator $set_seqSt_update_neighbSt",{
  
  # Initialize singleStructureGenerator instance
  single_obj <- singleStructureGenerator$new("U",10)
  
  # Expect error when either index or newSt arguments are not given
  expect_error(single_obj$set_seqSt_update_neighbSt(newSt = 2),
               info = "method fails to throw error when 'index' argument is not given")
  expect_error(single_obj$set_seqSt_update_neighbSt(index = 2),
               info = "method fails to throw error when 'newSt' argument is not given")
  
  # Expect error when index is not one numerical value or not a value without decimals
  expect_error(single_obj$set_seqSt_update_neighbSt(index = c(1,10), newSt = 2),
               info = "method fails to throw error when 'index' argument is not one numerical value")
  expect_error(single_obj$set_seqSt_update_neighbSt(index = 1.5, newSt = 2),
               info = "method fails to throw error when 'index' argument has decimals")
  
  # Expect error when newSt is not 1, 2 or 3 (for unmethylated, partially-methylated or methylated)
  expect_error(single_obj$set_seqSt_update_neighbSt(index = 1, newSt = 4),
               info = "method fails to throw error when 'newSt' is not 1, 2 or 3")
  
  # Expect error when index given does not exist within the singleStructure instance length
  expect_error(single_obj$set_seqSt_update_neighbSt(index = 11, newSt = 2),
               info = "method fails to throw error when 'index' is not within sequence length")
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(single_obj$set_seqSt_update_neighbSt(index = 2, newSt = 2),
              info = "whith testing = FALSE method generates output")
  
  # Test sequence state is correctly changed
  output <- single_obj$set_seqSt_update_neighbSt(index = 4, newSt = 3, testing = TRUE)
  expect_equal(output$seq[4], 3,
               info = "method fails to assign new methylaton state to the given index")
  
  # Test neighbSt is correctly updated
  mapNeighbSt_matrix = matrix(c(1L:9L), byrow = TRUE, nrow = 3)
  exp3rdsite_neighbSt <- mapNeighbSt_matrix[output$seq[2], output$seq[4]]
  exp5thsite_neighbSt <- mapNeighbSt_matrix[output$seq[4], output$seq[6]]
  expect_equal(output$neighbSt[3], exp3rdsite_neighbSt,
               info = "method fails to update neighbSt to left neighbor")
  expect_equal(output$neighbSt[5], exp5thsite_neighbSt,
               info = "method fails to update neighbSt to right neighbor")
})
