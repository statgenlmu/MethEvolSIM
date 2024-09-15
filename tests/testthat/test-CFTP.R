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

test_that("singleStructureGenerator $get_seqSt_leftneighb() and $get_seqSt_rightneighb()",{
  
  # Initialize singleStructureGenerator instance
  single_obj <- singleStructureGenerator$new("U",100)
  
  # Expect error when function is called without giving a value for the index 
  expect_error(single_obj$get_seqSt_leftneighb(),
               info = "method fails to throw error when no value for 'index' argument is given")
  expect_error(single_obj$get_seqSt_rightneighb(),
               info = "method fails to throw error when no value for 'index' argument is given")
  
  # Expect error when function is called with a non-valid index value
  expect_error(single_obj$get_seqSt_leftneighb(index = "U"),
               info = "method fails to throw error when 'index' argument is not one integer value")
  expect_error(single_obj$get_seqSt_rightneighb(index = "U"),
               info = "method fails to throw error when 'index' argument is not one integer value")
  expect_error(single_obj$get_seqSt_leftneighb(index = c(1,2)),
               info = "method fails to throw error when 'index' argument is not one integer value")
  expect_error(single_obj$get_seqSt_rightneighb(index = c(1,2)),
               info = "method fails to throw error when 'index' argument is not one integer value")
  expect_error(single_obj$get_seqSt_leftneighb(index = 1.5),
               info = "method fails to throw error when 'index' argument is not one integer value")
  expect_error(single_obj$get_seqSt_rightneighb(index = 1.5),
               info = "method fails to throw error when 'index' argument is not one integer value")
  
  # Expect error when function is called with an index value not within sequence length
  expect_error(single_obj$get_seqSt_leftneighb(index = -1),
               info = "method fails to throw error when 'index' argument is not within sequence length")
  expect_error(single_obj$get_seqSt_rightneighb(index = 200),
               info = "method fails to throw error when 'index' argument is not within sequence length")
  
  # Test that the sequence state for each site's left and right neighbor is correct
  for (s in 2:10){
    expect_equal(single_obj$get_seqSt_leftneighb(s), single_obj$get_seq()[s-1],
                 info = "sequence methylation state obtained with get_leftneighbSt is not correct")
  }
  for (s in 1:9){
    expect_equal(single_obj$get_seqSt_rightneighb(s), single_obj$get_seq()[s+1],
                 info = "neighbSt obtained with get_rightneighbSt is not correct")
  }
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

test_that("combiStructureGenerator $cftp_event_generator()", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  
  # Expect error when steps argument is not given
  expect_error(combi_obj$cftp_event_generator(),
               info = "method fails to throw error when 'steps' argument is not given")
  
  # Expect error when steps argument is not an non-decimal numerical value of min value 1
  expect_error(combi_obj$cftp_event_generator(steps = "5"),
               info = "method fails to throw error when 'steps' argument is non-numeric")
  expect_error(combi_obj$cftp_event_generator(steps = c(1,5)),
               info = "method fails to throw error when 'steps' argument is not of length 1")
  expect_error(combi_obj$cftp_event_generator(steps = 0),
               info = "method fails to throw error when 'steps' argument is < 1")
  expect_error(combi_obj$cftp_event_generator(steps = 100.5),
               info = "method fails to throw error when 'steps' argument has decimal numbers")
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(combi_obj$cftp_event_generator(steps = 1000),
              info = "whith testing = FALSE method generates output")
  
  # Test proposed event, singleStr and site, and threshold are generated for each CFTP step
  combi_obj <- combiStructureGenerator$new(infoStr)
  
  # 1: Test length of events after calling the method once and twice
  example_steps <- 1000
  output <- combi_obj$cftp_event_generator(steps = example_steps, testing = TRUE)
  expect_equal(length(output$CFTP_chosen_singleStr), example_steps,
               info = "length of CFTP_chosen_singleStr not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_chosen_site), example_steps,
               info = "length of CFTP_chosen_site not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_event), example_steps,
               info = "length of CFTP_event not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_random), example_steps,
               info = "length of CFTP_random not equal to number of steps after first event generation")
  
  example_steps_2 <- 50
  output <- combi_obj$cftp_event_generator(steps = example_steps_2, testing = TRUE)
  total_steps <- example_steps + example_steps_2
  expect_equal(length(output$CFTP_chosen_singleStr), total_steps,
               info = "length of CFTP_chosen_singleStr not equal to number of steps after second event generation")
  expect_equal(length(output$CFTP_chosen_site), total_steps,
               info = "length of CFTP_chosen_site not equal to number of steps after second event generation")
  expect_equal(length(output$CFTP_event), total_steps,
               info = "length of CFTP_event not equal to number of steps after first second generation")
  expect_equal(length(output$CFTP_random), total_steps,
               info = "length of CFTP_random not equal to number of steps after first second generation")
  
  # 2: Test content of events as expected
  expect_true(all(output$CFTP_chosen_singleStr %in% c(1,2,3)),
              info = "CFTP_chosen_singleStr contains singleStr index not in combiStr")
  expect_true(all(output$CFTP_chosen_site %in% 1:10),
              info = "CFTP_chosen_site contains site indexes not in combiStr")
  expect_true(all(output$CFTP_event %in% 1:5),
              info = "CFTP_event contains events not in 1:5")
  expect_true(all(output$CFTP_random >= 0 & output$CFTP_random <= 1),
              info = "CFTP_random contains threshold not between 0 and 1")
  
  #3: Test case with one singleStr
  infoStr <- data.frame(n = c(10),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  example_steps <- 100
  output <- combi_obj$cftp_event_generator(steps = example_steps, testing = TRUE)
  
  expect_true(all(output$CFTP_chosen_singleStr == 1),
              info = "CFTP_chosen_singleStr contains singleStr index not in combiStr")
  expect_true(all(output$CFTP_chosen_site %in% 1:10),
              info = "CFTP_chosen_site contains site indexes not in combiStr")

  
  #3: Test case with one singleStr with one position
  infoStr <- data.frame(n = c(1),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  example_steps <- 100
  output <- combi_obj$cftp_event_generator(steps = example_steps, testing = TRUE)
  
  expect_true(all(output$CFTP_chosen_singleStr == 1),
              info = "CFTP_chosen_singleStr contains singleStr index not in combiStr")
  expect_true(all(output$CFTP_chosen_site ==1),
              info = "CFTP_chosen_site contains site indexes not in combiStr")

})

test_that("combiStructureGenerator $cftp_apply_events()", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  
  # Expect error when calling method before using $cftp_event_generator() to generate events
  expect_error(combi_obj$cftp_apply_events(),
               info = "method fails to throw error when called before generating CFTP events")
  
  # Test output
  combi_obj$cftp_event_generator(steps = 100)
  output <- combi_obj$cftp_apply_events(testing = TRUE)
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(combi_obj$cftp_apply_events(),
              info = "whith testing = FALSE method generates output")
})


