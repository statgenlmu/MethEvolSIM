test_that("cftpStepGenerator_new $new", {
  
  cftp <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  
  expect_equal(cftp$singleStr_number, 3,
               info = "Incorrect number of singleStr instances")
  expect_equal(cftp$singleStr_siteNumber, c(10, 8, 15),
               info = "Incorrect number of sites in singleStr instances")
  expect_equal(cftp$CFTP_highest_rate,  0.8,
               info = "Incorrect value of CFTP highest rate")
  expect_equal(cftp$number_steps, 0,
               info = "Incorrect number of steps")
  expect_true(is.null(cftp$steps_perVector),
              info = "Not null value of steps per vector at initialization")
  
  expect_true(is.list(cftp$CFTP_chosen_singleStr),
              info = "Attribute $CFTP_chosen_singleStr not initialized as list")
  expect_equal(length(cftp$CFTP_chosen_singleStr), 0,
               info = "Attribute $CFTP_chosen_singleStr not initialized with length 0")
  
  expect_true(is.list(cftp$CFTP_chosen_site),
              info = "Attribute $CFTP_chosen_site not initialized as list")
  expect_equal(length(cftp$CFTP_chosen_site), 0,
               info = "Attribute $CFTP_chosen_site not initialized with length 0")
  
  expect_true(is.list(cftp$CFTP_event),
              info = "Attribute $CFTP_event not initialized as list")
  expect_equal(length(cftp$CFTP_event), 0,
               info = "Attribute $CFTP_event not initialized with length 0")
  
  expect_true(is.list(cftp$CFTP_random),
              info = "Attribute $CFTP_random not initialized as list")
  expect_equal(length(cftp$CFTP_random), 0,
               info = "Attribute $CFTP_random not initialized with length 0")
  
})

test_that("cftpStepGenerator_new $generate_events", {
  
  cftp <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  
  # Expect error when steps argument is not an non-decimal numerical value of min value 1
  expect_error(cftp$generate_events(steps = "5"),
               info = "method fails to throw error when 'steps' argument is non-numeric")
  expect_error(cftp$generate_events(steps = c(1,5)),
               info = "method fails to throw error when 'steps' argument is not of length 1")
  expect_error(cftp$generate_events(steps = 0),
               info = "method fails to throw error when 'steps' argument is < 1")
  expect_error(cftp$generate_events(steps = 100.5),
               info = "method fails to throw error when 'steps' argument has decimal numbers")
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(cftp$generate_events(steps = 1000),
              info = "whith testing = FALSE method generates output")
  
  ##TODO:HERE add testing of old_steps == $number_steps and new_steps and current_listIndex and self$steps_perVector
  # 1: Test length of events after calling the method once and twice
  example_steps <- 1000
  cftp <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  output <- cftp$generate_events(steps = example_steps, testing = TRUE)
  expect_equal(length(output$CFTP_chosen_singleStr), example_steps,
               info = "length of CFTP_chosen_singleStr not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_chosen_site), example_steps,
               info = "length of CFTP_chosen_site not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_event), example_steps,
               info = "length of CFTP_event not equal to number of steps after first event generation")
  expect_equal(length(output$CFTP_random), example_steps,
               info = "length of CFTP_random not equal to number of steps after first event generation")
  # Expect error if called again with smaller number of steps
  expect_error(cftp$generate_events(steps = 50),
               info = "method fails to stop when given steps have already been generated")
  
  total_steps <- 2000
  output <- cftp$generate_events(steps = total_steps, testing = TRUE)
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
  expect_true(all(output$CFTP_chosen_site %in% 1:15),
              info = "CFTP_chosen_site contains site indexes not in combiStr")
  expect_true(all(output$CFTP_event %in% 1:5),
              info = "CFTP_event contains events not in 1:5")
  expect_true(all(output$CFTP_random >= 0 & output$CFTP_random <= 1),
              info = "CFTP_random contains threshold not between 0 and 1")
})