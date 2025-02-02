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



test_that("cftpStepGenerator $generate_events according to combi instance", {
  
  # test 1
  test <- 1
  infoStr <- data.frame(n = c(25, 78, 1, 16, 88, 40),
                        globalState = c("M", "U", "M", "M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  cftp <- cftpStepGenerator$new(singleStr_number = c$get_singleStr_number(),
                                singleStr_siteNumber = c$get_singleStr_siteNumber(), 
                                CFTP_highest_rate = c$get_highest_rate())
  cftp$generate_events()

  expect_true(all(cftp$CFTP_chosen_singleStr %in% 1:6),
              info = paste("Samples singleStr indices not in combistrucutre instance in test:", test))
  
  for (str in 1:6){
    singleStr_indices <- which(cftp$CFTP_chosen_singleStr == str)
    expect_true(all(cftp$CFTP_chosen_site[singleStr_indices] %in% 1:infoStr$n[str]),
                info = paste("Samples site indice not within singleStr length in test:", test, "singleStr index:", str))
  }
  
    
  # test 2: one singleStr
  test <- 2
  infoStr <- data.frame(n = c(10),
                        globalState = c("M"))
  c <- combiStructureGenerator$new(infoStr)
  cftp <- cftpStepGenerator$new(singleStr_number = c$get_singleStr_number(),
                                singleStr_siteNumber = c$get_singleStr_siteNumber(), 
                                CFTP_highest_rate = c$get_highest_rate())
  cftp$generate_events()
  
  expect_true(all(cftp$CFTP_chosen_singleStr == 1),
              info = paste("Samples singleStr indices not in combistrucutre instance in test:", test))
  expect_true(all(cftp$CFTP_chosen_site %in% 1:10),
              info = paste("Samples site indices not in combistrucutre instance in test:", test))
  
  # test 3: one singleStr with one position
  test <- 3
  infoStr <- data.frame(n = c(1),
                        globalState = c("M"))
  c <- combiStructureGenerator$new(infoStr)
  cftp <- cftpStepGenerator$new(singleStr_number = c$get_singleStr_number(),
                                singleStr_siteNumber = c$get_singleStr_siteNumber(), 
                                CFTP_highest_rate = c$get_highest_rate())
  cftp$generate_events()
  
  expect_true(all(cftp$CFTP_chosen_singleStr == 1),
              info = paste("Samples singleStr indices not in combistrucutre instance in test:", test))
  expect_true(all(cftp$CFTP_chosen_site == 1),
              info = paste("Samples site indices not in combistrucutre instance in test:", test))
  
  
})

test_that("combiStructureGenerator $set_CFTP_info", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  # Expect error if given object is not of class cftpStepGenerator
  expect_error(c$set_CFTP_info(1),
               info = "method fails to throw an error when given object is not instance of class cftpStepGenerator")
  
  # Initialize cftpStepGenerator instance
  cftp <- cftpStepGenerator$new(singleStr_number = c$get_singleStr_number(),
                                singleStr_siteNumber = c$get_singleStr_siteNumber(), 
                                CFTP_highest_rate = c$get_highest_rate())
  
  c$set_CFTP_info(cftp)
  expect_true(class(c$get_CFTP_info())[1] == "cftpStepGenerator")
  expect_equal(length(c$get_CFTP_info()$CFTP_event), 0,
               info = "Assigned cftp instance without steps has length of CFTP_event non 0")
  
  cftp$generate_events(100)
  expect_equal(length(c$get_CFTP_info()$CFTP_event), 100,
               info = "Assigned cftp instance with 100 steps has length of CFTP_event non 100")
})


test_that("combiStructureGenerator $cftp_apply_events()", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  
  # Expect error when calling method before generating events
  expect_error(c$cftp_apply_events(),
               info = "method fails to throw error when called before generating CFTP events")
  
  # Test output
  
  # Set correspondence of event number and applied case
  encoding_case <- data.frame(n = 1:5,
                              case = c("SSEi_1", "SSEi_2", "SSEi_3", "SSEc_left", "SSEc_right"))
  
  # Initialize cftpStepGenerator instance and assign it to combiStructureGenerator instance
  cftp <- cftpStepGenerator$new(singleStr_number = c$get_singleStr_number(),
                                singleStr_siteNumber = c$get_singleStr_siteNumber(), 
                                CFTP_highest_rate = c$get_highest_rate())
  c$set_CFTP_info(cftp)
  cftp$generate_events(steps = 100)
  output <- c$cftp_apply_events(testing = TRUE)
  
  accepted_indeces <- which(output$event_acceptance == TRUE)
  if (length(accepted_indeces) > 0){
    expect_true(all(output$r_jk[accepted_indeces]/output$r_m > output$CFTP_random[accepted_indeces]),
                info = "Not all accepted events fulfill relative rate higher than sampled threshold")
    expect_false(any(is.na(output$CFTP_event[accepted_indeces])),
                 info = "Not all cases of accepted events return non-NA value in testing output applied_event")
    for(e in 1:length(accepted_indeces)){
      expect_equal(encoding_case[output$CFTP_event[accepted_indeces][e], "case"], output$applied_event[accepted_indeces][e],
                   info = "Error in correspondence between CFTP_event encoding and applied event")
    }
    
  }
  
  non_accepted_indeces <- which(output$event_acceptance == FALSE)
  if (length(non_accepted_indeces) > 0){
    expect_true(all(output$CFTP_event[non_accepted_indeces][is.na(output$r_jk[non_accepted_indeces])] %in% c(1, 2, 3)),
                info = "Not all cases in which a rate was not sampled because of SSEi newSt and oldSt being equal correspond to SSEi events")
    expect_true(all(output$r_jk[non_accepted_indeces][!is.na(output$r_jk[non_accepted_indeces])]/output$r_m <= output$CFTP_random[non_accepted_indeces][!is.na(output$r_jk[non_accepted_indeces])]),
                info = "Not all non_accepted_events with rate fulfill relative rate smaller or equal than sampled threshold")
    expect_true(all(is.na(output$applied_event[non_accepted_indeces])),
                info = "Not all cases of non-accepted events return NA in testing output applied_event")
    
  }
  
  # Expect NULL output when arguments are correct but testing is (as default) FALSE
  expect_null(c$cftp_apply_events(),
              info = "whith testing = FALSE method generates output")
})

##TODO: delete this if moved to the multiRegion script. or move to the top if this is kept as separate script
# Function to access private variables and functions
get_private <- function(x) {
  x[['.__enclos_env__']]$private
}

test_that("combiStructureGenerator $cftp", {
  
  # Initialize combiStructureGenerator instance
  infoStr <- data.frame(n = c(10, 10, 10),
                        globalState = c("M", "U", "M"))
  combi_obj <- combiStructureGenerator$new(infoStr)
  
  # Test total CFTP step number is minimum the given number
  test_steps <- 100
  output<- combi_obj$cftp(steps = test_steps, testing = TRUE)
  expect_true(output$total_steps >= test_steps,
              info = "method testing output total_steps is smaller than given one")
  
  # Extract sequences of combi instances to test they are equal
  m_seq <- c(output$combi_m$get_singleStr(1)$get_seq(), output$combi_m$get_singleStr(2)$get_seq(), output$combi_m$get_singleStr(3)$get_seq())
  self_seq <- c(output$self$get_singleStr(1)$get_seq(), output$self$get_singleStr(2)$get_seq(), output$self$get_singleStr(3)$get_seq())
  expect_true(all(m_seq == self_seq),
              info = "method testing output with different $seq in combi_u and self")
  
  # Expect method to return nothing when testing is (as default) FALSE
  combi_obj <- combiStructureGenerator$new(infoStr)
  expect_null(combi_obj$cftp(steps = test_steps),
               info = "method returns something when testing is FALSE")
  
  # Expect objects with correct combiStructure ID
  # Initiate an instance with a reset shared counter
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  # Reset shared counter
  c$reset_sharedCounter()
  c <- combiStructureGenerator$new(infoStr)
  output <- c$cftp(testing = TRUE)
  expect_equal(output$self$get_id(), 1, 
               info = "self ID not equal to one after resetting counter")
  expect_equal(output$combi_m$get_id(), output$counter+1,
               info = "combi_m ID not equal to 2* the counter of cftp cycles before convergence")
  
  for (i in 1:3){
    expect_equal(get_private(output$self$get_singleStr(i))$my_combiStructure$get_id(), 1,
                 info = paste("self singleStr", i, "does not point to correct my_combiStructure according to self ID"))
    expect_equal(get_private(output$combi_m$get_singleStr(i))$my_combiStructure$get_id(), output$counter+1,
                 info = paste("combi_m singleStr", i, "does not pont to correct my_combiStructure according to combi_m ID"))
  }
  
  # Expect number of cftp steps corresponding to default_steps*2^(counter-1)
  default_steps <- 10000
  calculate_total_steps <- function(counter) {
    10000 * 2^(counter - 1)
  }
  expect_equal(output$total_steps, calculate_total_steps(output$counter),
               info = "incorrect number of steps generated for given number of cftp cycles before convergence")
  
  
  # Test case with several singleStr with different lengths
  infoStr <- data.frame(n = c(10, 35, 1),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  output <- c$cftp(testing = TRUE)
  
  expect_true(all(sort(unique(output$CFTP_chosen_site))==1:35),
              info = "not all possible sites are chosen when singleStructures have different lengths")
  
  # Expect clones to have the same CFTP_info
  expect_true(identical(output$self$get_CFTP_info(), output$combi_m$get_CFTP_info()),
              info = "Clones have non-identical CFTP_info")
})

test_that("treeMultiRegionSimulator initialize with cftp = TRUE", {
  # Initiate a combiStructure instance to reset shared counter
  infoStr <- data.frame(n = c(13, 13, 13),
                        globalState = c("M", "U", "M"))
  c <- combiStructureGenerator$new(infoStr)
  # Reset shared counter
  c$reset_sharedCounter()
  
  # Initialize treeMultiRegionSimulator instance
  infoStr <- data.frame(n = c(20, 10, 20),
                        globalState = c("M", "U", "M"))
  message <- capture.output(t <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1);", CFTP = TRUE, testing = TRUE), type = "message")
  
  # Extract sequence info
  root_before_cftp_seq <- t$testing_output$seq_before_cftp
  root_after_cftp_seq <- c()
  self_seq <- c()
  for (i in 1:3){
    root_after_cftp_seq <- c(root_after_cftp_seq, t$Branch[[1]]$get_singleStr(i)$get_seq())
    self_seq <- c(self_seq, t$testing_output$cftp_output$self$get_singleStr(i)$get_seq())
  }
  
  expect_true(all(unique(sort(t$testing_output$cftp_output$CFTP_chosen_site)) == 1:20),
              info = "Not all sites chosen by $cftp")
  # Expect sequence changing after cftp and root sequence to be as in cftp output
  # Not always
  #expect_false(all(root_after_cftp_seq == root_before_cftp_seq),
  #             info = "$seq doesnt change after cftp")
  expect_true(all(self_seq == root_after_cftp_seq),
              info = "Root $seq doesnt correspond to $seq after cftp")
  
  # Expect ID in root data to be 1
  expect_equal(t$Branch[[1]]$get_id(), 1,
               info = "root ID does not stay as initial ID")
  
  # Expect new ID in root data to correspond to cftp's combi_u id
  expect_equal(t$testing_output$cftp_output$combi_m$get_id(), t$testing_output$cftp_output$counter +1,
               info = "$cftp combi_m ID does not correspond to the counter")
  for (i in 1:3){
    expect_equal(get_private(t$Branch[[1]]$get_singleStr(i))$my_combiStructure$get_id(), 1,
                 info = paste("root singleStr", i, "does not point to correct my_combiStructure according to root ID"))
    expect_equal(get_private(t$testing_output$cftp_output$combi_m$get_singleStr(i))$my_combiStructure$get_id(), t$testing_output$cftp_output$counter +1,
                 info = paste("combi_m singleStr", i, "does not point to correct my_combiStructure according to counter"))
  }
  
  # Test ID change consistent with testing = FALSE
  c$reset_sharedCounter()
  message <- capture.output(t <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1);", CFTP = TRUE), type = "message")
  
  expect_true(t$Branch[[1]]$get_id() == 1,
              info = "id with testing = FALSE remains consistent to testing = TRUE")
  
  
  # Expect informative message when CFTP = TRUE and no message about CFTP when is, as default, FALSE
  message <- capture.output(t <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1);", CFTP = TRUE), type = "message")
  
  expect_true(message[2] == "Calling CFTP algorithm for data at root before letting it evolve along given tree.",
              info = "Class fails to inform when CFTP is set to TRUE")
  
  message <- capture.output(t <- treeMultiRegionSimulator$new(infoStr = infoStr, tree = "(a:1);"), type = "message")
  
  expect_true(is.na(message[2]),
              info = "Class contains informative message when CFTP is set to FALSE")
})

##TODO: comments to week_6
##TODO: prepare script to run a case with a seed with current branch and another with previous branch
##TODO: prepare a script to compare the output with both versions, which should be identical
##TODO: start in server previously crushed runs. If they dont crush anymore:
##TODO: clean old tests
##TODO: Update documentation
