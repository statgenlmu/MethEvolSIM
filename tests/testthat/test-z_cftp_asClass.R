test_that("combiStructureGenerator $get_highest_rate", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  output <- c$get_highest_rate()
  
  expect_true(all(is.numeric(output), length(output) == 1, output > 0))
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

test_that("cftpStepGenerator $new", {
  
  ## TODO: Add test for checking that when one argument is not given, the method throws an error
  
  cftp <- cftpStepGenerator$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
  
  expect_equal(cftp$singleStr_number, 3,
               info = "Incorrect number of singleStr instances")
  expect_equal(cftp$singleStr_siteNumber, c(10, 8, 15),
               info = "Incorrect number of sites in singleStr instances")
  expect_equal(cftp$CFTP_highest_rate,  0.8,
              info = "Incorrect value of CFTP highest rate")
  
})

test_that("combiStructureGenerator $cftp initialization of cftpStepGenerator instance", {
  
  # test 1
  test <- 1
  infoStr <- data.frame(n = c(10, 8, 15),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  expect_true(is.null(c$get_CFTP_info()),
              info = paste("Not null info in combi instance before calling $cftp() in test:", test))
  
  c$cftp()
  
  expect_equal(class(c$get_CFTP_info())[1], "cftpStepGenerator")
  expect_equal(c$get_CFTP_info()$singleStr_number, length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(c$get_CFTP_info()$singleStr_siteNumber, infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
  clone <- c$copy()
  expect_equal(clone$get_CFTP_info()$singleStr_number, length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(clone$get_CFTP_info()$singleStr_siteNumber, infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
  
  # test 2
  infoStr <- data.frame(n = c(25, 78, 1, 16, 88, 40),
                        globalState = c("M", "U", "M", "M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  expect_true(is.null(c$get_CFTP_info()),
              info = paste("Not null info in combi instance before calling $cftp() in test:", test))
  
  c$cftp()
  
  expect_equal(class(c$get_CFTP_info())[1], "cftpStepGenerator")
  expect_equal(c$get_CFTP_info()$singleStr_number, length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(c$get_CFTP_info()$singleStr_siteNumber, infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
  clone <- c$copy()
  expect_equal(clone$get_CFTP_info()$singleStr_number, length(infoStr$n),
               info = paste("Incorrect number of singleStr instances in test:", test))
  expect_equal(clone$get_CFTP_info()$singleStr_siteNumber, infoStr$n,
               info = paste("Incorrect number of sites in singleStr instances in test:", test))
  
})

test_that("cftpStepGenerator $generate_events", {
  
  cftp <- cftpStepGenerator$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)

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
  
  # 1: Test length of events after calling the method once and twice
  example_steps <- 1000
  cftp <- cftpStepGenerator$new(singleStr_number = 3, singleStr_siteNumber = c(10, 8, 15), CFTP_highest_rate = 0.8)
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

test_that("cftpStepGenerator $generate_events according to combi instance", {
  
  # test 1
  infoStr <- data.frame(n = c(25, 78, 1, 16, 88, 40),
                        globalState = c("M", "U", "M", "M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)

    
  # test 2: one singleStr
  infoStr <- data.frame(n = c(10),
                        globalState = c("M"))
  c <- combiStructureGenerator$new(infoStr)
  
  # test 3: one singleStr with one position
  infoStr <- data.frame(n = c(1),
                        globalState = c("M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  
  #TODO: Complete using c$cftp when method is ready by modifying last two cases in test combiStructureGenerator $cftp_event_generator()
})

##TODO: next update cftp_apply_events()
