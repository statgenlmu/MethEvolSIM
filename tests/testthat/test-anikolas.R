library(devtools)

#Recall: Xseq means everything from the left of the recombination point and Yseq means everything to the right

get_private <- function(x){
  x[['.__enclos_env__']]$private
}

#testing get_SrMatrix properties
#when a property is violated a warning gets thrown by get_SrMatrix
test_that("get_TransMatrix property test",{
test_num <- 100
#iterating through (most likely) different test cases
  for (x in 1:test_num){
    #creating random vectors with length 3 summing to 1
    random_vector1 <- runif(3) / sum(runif(3))
    random_vector1 <- random_vector1 /sum(random_vector1)
    random_vector2 <- runif(3) / sum(runif(3))
    random_vector2 <- random_vector2 /sum(random_vector2)
    
    old_freqs <- random_vector1
    new_freqs <- random_vector2
    
    if(all.equal(old_freqs, new_freqs)!=TRUE){
      #creating a new ssg object to test get_SrMatrix
      ssg <- singleStructureGenerator$new(globalState = "U", n = x,testing=FALSE)
      #testing is implemented in get_SrMatrix so we expect no warnings
      expect_no_warning(ssg$get_TransMatrix(old_eqFreqs = old_freqs, new_eqFreqs = new_freqs, testing=TRUE))
    }
    else{
      test_num <- test_num + 1
    }
  }
})

#testing RE_within properties with random vectors
#when a property is violated a warning gets thrown by RE_within
test_that("RE_within property test",{
  test_num <- 100
  for(i in 5:test_num){
    vec <- sample(c(1, 2, 3), i-3, replace = TRUE)
    #creating random vectors summing to 1
    random_vector3 <- runif(3) / sum(runif(3))
    random_vector3 <- random_vector3 /sum(random_vector3)
    random_vector4 <- runif(3) / sum(runif(3))
    random_vector4 <- random_vector4 /sum(random_vector4)
    if(all.equal(random_vector3,random_vector4)!= TRUE){
      #test if RE_within throws any warnings
      ssg <- singleStructureGenerator$new(globalState = "U", n = i,testing=FALSE)
      expect_no_warning(ssg$RE_within(Y_seq = vec,Y_eqFreqs = random_vector3,testing=TRUE))
      
      #test whether the chosen frequencies in RE_within are the same as the freqs set in the ssg
      expect_equal(ssg$RE_within(Y_seq = vec,Y_eqFreqs = random_vector3,testing=TRUE)$chosenEqFreqs, ssg$get_eqFreqs())
    }
    else{
      test_num <- test_num+1
    }
  }
})

#test whether the eqfrequencies calculated in RE_within make sense by recalculating them afterwards for the FreqsX part
#this only works when the size of the sequence is large. with smaller sizes the test may fail
test_that("RE_within frequency test",{
  #different testing cases
  testing_parameters <- list(list(Y_seq = c(rep(1, 1500)), yEqFreqs = c(0.9,0.09,0.01), globalState = "U",n=3000),
                             list(Y_seq = c(rep(2,500)), yEqFreqs = c(0.05,0.9,0.05), globalState = "U",n=2000),
                             list(Y_seq = c(rep(3,100)), yEqFreqs = c(0,0.01,0.99), globalState = "U",n=1000),
                             list(Y_seq = c(rep(2,200)), yEqFreqs = c(0,0.99,0.01), globalState = "M",n=400)
                             )
  #iterating through each test case and extracting the information to call RE_within with testing=TRUE
  for(i in 1:length(testing_parameters)){
    Y_seq_test <- testing_parameters[[i]]$Y_seq
    yEqFreqs_test <- testing_parameters[[i]]$yEqFreqs
    globalState_test <- testing_parameters[[i]]$globalState
    n_test  <- testing_parameters[[i]]$n
    ssg <- singleStructureGenerator$new(globalState=globalState_test,n=n_test)
    xEqFreqs_test <- ssg$get_eqFreqs()
    out <- ssg$RE_within(Y_seq= Y_seq_test, Y_eqFreqs = yEqFreqs_test,testing = TRUE)
    
    #reminder: X_Seq is the left part of the sequence relative to the recombination point
    #and Y_Seq is the right part relative to the recombination point
    
    #when the FreqsX were chosen we calculate the expected frequencies for the Y_Seq part
    if (out$id =="FreqsX"){
      xSeq <- out$chosenSeq
      len_x <- length(xSeq)
      len_y <- length(ssg$get_seq())-len_x
      
      pi_u_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==1))/len_y
      pi_p_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==2))/len_y
      pi_m_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==3))/len_y
      expect_equal(c(pi_u_x,pi_p_x,pi_m_x),xEqFreqs_test,tolerance=0.2)
    }
    #when the FreqsY were chosen we calculate the expected frequencies for the xseq part
    if (out$id == "FreqsY"){
      len <- length(Y_seq_test)
      pi_u <- length(which(ssg$get_seq()[1:len]==1))/len
      pi_p <- length(which(ssg$get_seq()[1:len]==2))/len
      pi_m <- length(which(ssg$get_seq()[1:len]==3))/len
      expect_equal(c(pi_u,pi_p,pi_m),yEqFreqs_test,tolerance = 0.1)
    }
  }
})

#here i test whether the X and Y sequence that is not part of the recombining ssg stays the same after recombination as it should
test_that("RE_complete unchanged X/Y sequence test",{
  #testing cases
  testing_parameters <- list(list(infoStr_1 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),position = 30),
                             list(infoStr_1 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),infoStr_2 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),position = 3),
                             list(infoStr_1 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),position = 3),
                             list(infoStr_1 = data.frame(n = c(1000,10),globalState = c("U","U")),infoStr_2 = data.frame(n = c(1000,10),globalState = c("M","M")),position = 1004)
                            )
  #iterating through every test case
  for(i in 1:length(testing_parameters)){
    #testing for left part of ssgs to the recombining ssg
    infoStr_1_testing <- testing_parameters[[i]]$infoStr_1
    infoStr_2_testing <- testing_parameters[[i]]$infoStr_2
    position_testing <- testing_parameters[[i]]$position
    
    #create two combistructure objects and clone them to be able to compare old and new combistructure after recombinaion
    csg_1 <- combiStructureGenerator$new(infoStr = infoStr_1_testing)
    csg_1_copy <- csg_1$copy()
    csg_2 <- combiStructureGenerator$new(infoStr = infoStr_2_testing)
    csg_2_copy <- csg_2$copy()
    
    #recombination event on csg_1 and csg_2 at position position_testing
    out <- csg_1$RE_complete(combi_2 = csg_2,position = position_testing,testing=TRUE)
    ssg_position <- out$ssg_index
    ssg_within_position <- out$ssg_position_within
    
    #this is the sequence of ssgs before the recombination event if the Xseq was chosen
    seq_before_x <- c()
    #this is the sequence of ssgs after the recombination event if the Xseq was chosen
    seq_after_x <- c()

    eqFreqs_before_x <- c()
    eqFreqs_after_x <- c()
    
    #Iterating through every ssg index and appending the corresponding ssg to seq_before_x and seq_after_x.
    #We expect seq_before_x and seq_after_x to be equal
    for(i in 1:ssg_position){
      
      if(i == ssg_position){
        #the correctness of the recombining ssg is already tested in the RE_within test, we are only interested in the unaffected single structures in this test
        next
        #Interesting thought:
        #these two lines down below are used for appending the sequence in the recombining ssg but i would need to include additional
        #information about which sequence got chosen during the RE_within process and that was already tested above so i skipped it
        #update: it probably wouldnt make sense since the neighborstates get updated and thus the sequences differ from the earlier ones
        #seq_before_x <- append(seq_before_x, csg_1_copy$get_singleStr(i)$get_seq()[1:ssg_within_position-1])
        #seq_after_x <- append(seq_after_x,csg_1$get_singleStr(i)$get_seq()[1:ssg_within_position-1])
      }
      else{
        #appending the ssg sequence before and after recombination to a vector
        seq_before_x <- append(seq_before_x, csg_1_copy$get_singleStr(i)$get_seq())
        seq_after_x <- append(seq_after_x, csg_1$get_singleStr(i)$get_seq())
        
        eqFreqs_before_x <- append(eqFreqs_before_x, csg_1_copy$get_singleStr(i)$get_eqFreqs())
        eqFreqs_after_x <- append(eqFreqs_after_x, csg_1$get_singleStr(i)$get_eqFreqs())
      }
    }
    expect_equal(seq_before_x,seq_after_x,info="unaffected X sequence is not equal after recombination")
    expect_equal(eqFreqs_before_x,eqFreqs_after_x,info="unaffected X eqFreqs are not equal after recombination")
    
    #testing for the right part of ssgs after the recombining ssg
    #this is the sequence of ssgs before the recombination event if the Yseq was chosen
    seq_before_y <- c()
    #this is the sequence of ssgs after the recombination event if the Yseq was chosen
    seq_after_y <- c()
  
    eqFreqs_before_y <- c()
    eqFreqs_after_y <- c()
  
    for(i in ssg_position:csg_1$get_singleStr_number()){
      if(i==ssg_position){
        #this is already tested in the RE_within test
        next
        #Interesting thought:
        #these two lines down below are used for appending the sequence in the recombining ssg but i would need to include additional
        #information about which sequence got chosen during the RE_within process and that was already tested above so i skipped it
        #update: it probably wouldnt make sense since the neighborstates get updated and thus the sequences differ from the earlier ones
        #seq_before_y <- append(seq_before_y,csg_1$get_singleStr(i)$get_seq()[ssg_within_position:length(seq_before_y,csg_1$get_singleStr(i)$get_seq())])
        #seq_after_y <- append(seq_after_y, csg_2_copy$get_singleStr(i)$get_seq()[ssg_within_position:length(seq_before_y,csg_2$get_singleStr(i)$get_seq())])
      }
      else{
        seq_before_y <- append(seq_before_y, csg_1$get_singleStr(i)$get_seq())
        seq_after_y <- append(seq_after_y, csg_2_copy$get_singleStr(i)$get_seq())
      
        eqFreqs_before_y <- append(eqFreqs_before_y, csg_1$get_singleStr(i)$get_eqFreqs())
        eqFreqs_after_y <- append(eqFreqs_after_y, csg_2_copy$get_singleStr(i)$get_eqFreqs())
      }
    }
    expect_equal(seq_before_y,seq_after_y, info="unaffected Y sequence is not equal after recombination")
    expect_equal(eqFreqs_before_y,eqFreqs_after_y,info="unaffected Y eqFreqs are not equal after recombination")
  }
})

#testing whether the equilibrium frequencies get updated correctly after recombination events in a recombining combistructure object
test_that("RE_complete eqFreqs test",{
  #test cases
  testing_parameters <- list(list(infoStr_1 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),position = 30),
                             list(infoStr_1 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),infoStr_2 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),position = 3),
                             list(infoStr_1 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),position = 3),
                             list(infoStr_1 = data.frame(n = c(1000,10),globalState = c("U","U")),infoStr_2 = data.frame(n = c(1000,10),globalState = c("M","M")),position = 1004)
  )
  #iterating through the test cases
  for(i in 1:length(testing_parameters)){
    infoStr_1_testing <- testing_parameters[[i]]$infoStr_1
    infoStr_2_testing <- testing_parameters[[i]]$infoStr_2
    position_testing <- testing_parameters[[i]]$position
    
    csg_1 <- combiStructureGenerator$new(infoStr = infoStr_1_testing)
    csg_2 <- combiStructureGenerator$new(infoStr = infoStr_2_testing)
    
    #recombine the created combistructure objects
    out <- csg_1$RE_complete(combi_2 = csg_2,position = position_testing,testing=TRUE)
    ssg_position <- out$ssg_index
    ssg_eqFreqs <- out$recombEqFreqs
    
    #expect that the eqfreqs were correctly updated in the corresponding combistructre
    expect_equal(ssg_eqFreqs, csg_1$get_singleStr(ssg_position)$get_eqFreqs(),info="EqFreqs are not correctly applied to recombining ssg")
  }
})



















  
