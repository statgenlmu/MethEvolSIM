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

test_that("cftpStepGenerator_new $generate_events input control", {
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
})

test_that("generate_events correctly initializes steps", {
  obj <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  steps <- 10000
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(obj$generate_events(steps = steps),
              info = "whith testing = FALSE method generates output")
  
  obj <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  obj$generate_events(steps = steps)
  
  expect_equal(obj$steps_perVector, steps,
               info = "non-correct $steps_perVector after first generate_events call")
  expect_equal(obj$number_steps, steps,
               info = "non-correct number of steps after first generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr), 1,
               info = "non-correct number of list elements in $CFTP_chosen_singleStr after first generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr[[1]]), steps,
               info = "non-correct number of steps in $CFTP_chosen_singleStr after first generate_events call")
  expect_true(all(obj$CFTP_chosen_singleStr[[1]] %in% c(1,2,3)),
              info = "non-correct structure indices in $CFTP_chosen_singleStr after first generate_events call")
  expect_equal(length(obj$CFTP_chosen_site), 1,
               info = "non-correct number of list elements in $CFTP_chosen_site after first generate_events call")
  expect_equal(length(obj$CFTP_chosen_site[[1]]), steps,
               info = "non-correct number of steps in $CFTP_chosen_site after first generate_events call")
  expect_true(all(obj$CFTP_chosen_site[[1]] %in% 1:15),
              info = "non-correct site indices in $CFTP_chosen_site after first generate_events call")
  expect_equal(length(obj$CFTP_event), 1,
               info = "non-correct number of list elements in $CFTP_event after first generate_events call")
  expect_equal(length(obj$CFTP_event[[1]]), steps,
               info = "non-correct number of steps in $CFTP_event after first generate_events call")
  expect_true(all(obj$CFTP_event[[1]] %in% 1:5),
              info = "non-correct event encoding in $CFTP_event after first generate_events call")
  expect_equal(length(obj$CFTP_random), 1,
               info = "non-correct number of list elements in $CFTP_random after first generate_events call")
  expect_equal(length(obj$CFTP_random[[1]]), steps,
               info = "non-correct number of steps in $CFTP_random after first generate_events call")
  expect_true(all(obj$CFTP_random[[1]] >= 0 & obj$CFTP_random[[1]] <= 1),
              info = "non-correct values in $CFTP_random after first generate_events call")
})

test_that("generate_events prevents duplicate steps", {
  obj <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  obj$generate_events(steps = 10000)
  
  expect_error(obj$generate_events(steps = 5000), "The given number of steps has already been generated")
  expect_error(obj$generate_events(steps = 10000), "The given number of steps has already been generated")
})

test_that("generate_events handles step increments", {
  obj <- cftpStepGenerator_new$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  steps <- 10000
  obj$generate_events(steps = steps)
  
  # Increment steps as in second cftp iteration
  steps <- steps*2
  obj$generate_events(steps = steps)
  
  expect_equal(obj$steps_perVector, 10000,
               info = "non-correct $steps_perVector after second generate_events call")
  expect_equal(obj$number_steps, steps,
               info = "non-correct number of steps after second generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr), 2,
               info = "non-correct number of list elements in $CFTP_chosen_singleStr after second generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr[[2]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_singleStr after second generate_events call")
  expect_true(all(obj$CFTP_chosen_singleStr[[2]] %in% c(1,2,3)),
              info = "non-correct structure indices in $CFTP_chosen_singleStr after second generate_events call")
  expect_equal(length(obj$CFTP_chosen_site), 2,
               info = "non-correct number of list elements in $CFTP_chosen_site after second generate_events call")
  expect_equal(length(obj$CFTP_chosen_site[[2]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_site after second generate_events call")
  expect_true(all(obj$CFTP_chosen_site[[2]] %in% 1:15),
              info = "non-correct site indices in $CFTP_chosen_site after second generate_events call")
  expect_equal(length(obj$CFTP_event), 2,
               info = "non-correct number of list elements in $CFTP_event after second generate_events call")
  expect_equal(length(obj$CFTP_event[[2]]), 10000,
               info = "non-correct number of steps in $CFTP_event after second generate_events call")
  expect_true(all(obj$CFTP_event[[1]] %in% 1:5),
              info = "non-correct event encoding in $CFTP_event after second generate_events call")
  expect_equal(length(obj$CFTP_random), 2,
               info = "non-correct number of list elements in $CFTP_random after second generate_events call")
  expect_equal(length(obj$CFTP_random[[2]]), 10000,
               info = "non-correct number of steps in $CFTP_random after second generate_events call")
  expect_true(all(obj$CFTP_random[[2]] >= 0 & obj$CFTP_random[[2]] <= 1),
              info = "non-correct values in $CFTP_random after second generate_events call")
  
  # Increment steps as in third cftp iteration
  steps <- steps*2
  obj$generate_events(steps = steps)
  
  expect_equal(obj$steps_perVector, 10000,
               info = "non-correct $steps_perVector after third generate_events call")
  expect_equal(obj$number_steps, steps,
               info = "non-correct number of steps after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr), 4,
               info = "non-correct number of list elements in $CFTP_chosen_singleStr after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr[[3]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_singleStr after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_singleStr[[4]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_singleStr after third generate_events call")
  expect_true(all(obj$CFTP_chosen_singleStr[[3]] %in% c(1,2,3)),
              info = "non-correct structure indices in $CFTP_chosen_singleStr after third generate_events call")
  expect_true(all(obj$CFTP_chosen_singleStr[[4]] %in% c(1,2,3)),
              info = "non-correct structure indices in $CFTP_chosen_singleStr after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_site), 4,
               info = "non-correct number of list elements in $CFTP_chosen_site after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_site[[3]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_site after third generate_events call")
  expect_equal(length(obj$CFTP_chosen_site[[4]]), 10000,
               info = "non-correct number of steps in $CFTP_chosen_site after third generate_events call")
  expect_true(all(obj$CFTP_chosen_site[[3]] %in% 1:15),
              info = "non-correct site indices in $CFTP_chosen_site after third generate_events call")
  expect_true(all(obj$CFTP_chosen_site[[4]] %in% 1:15),
              info = "non-correct site indices in $CFTP_chosen_site after third generate_events call")
  expect_equal(length(obj$CFTP_event), 4,
               info = "non-correct number of list elements in $CFTP_event after third generate_events call")
  expect_equal(length(obj$CFTP_event[[3]]), 10000,
               info = "non-correct number of steps in $CFTP_event after third generate_events call")
  expect_equal(length(obj$CFTP_event[[4]]), 10000,
               info = "non-correct number of steps in $CFTP_event after third generate_events call")
  expect_true(all(obj$CFTP_event[[3]] %in% 1:5),
              info = "non-correct event encoding in $CFTP_event after third generate_events call")
  expect_true(all(obj$CFTP_event[[4]] %in% 1:5),
              info = "non-correct event encoding in $CFTP_event after third generate_events call")
  expect_equal(length(obj$CFTP_random), 4,
               info = "non-correct number of list elements in $CFTP_random after third generate_events call")
  expect_equal(length(obj$CFTP_random[[3]]), 10000,
               info = "non-correct number of steps in $CFTP_random after third generate_events call")
  expect_equal(length(obj$CFTP_random[[4]]), 10000,
               info = "non-correct number of steps in $CFTP_random after third generate_events call")
  expect_true(all(obj$CFTP_random[[3]] >= 0 & obj$CFTP_random[[3]] <= 1),
              info = "non-correct values in $CFTP_random after third generate_events call")
  expect_true(all(obj$CFTP_random[[4]] >= 0 & obj$CFTP_random[[4]] <= 1),
              info = "non-correct values in $CFTP_random after third generate_events call")
})


