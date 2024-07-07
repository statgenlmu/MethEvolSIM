library(devtools)
#testing get_SrMatrix
test_that("get_TransMatrix test",{
test_num <- 100
  for (x in 1:test_num){
    #creating random vectors with length 3 summing to 1
    random_vector1 <- runif(3) / sum(runif(3))
    random_vector1 <- random_vector1 /sum(random_vector1)
    random_vector2 <- runif(3) / sum(runif(3))
    random_vector2 <- random_vector2 /sum(random_vector2)
    
    old_freqs <- random_vector1
    new_freqs <- random_vector2
    
    if(all.equal(old_freqs, new_freqs)!=TRUE){
<<<<<<< HEAD
      #creating a new ssg object to test get_SrMatrix
      ssg <- singleStructureGenerator$new(globalState = "U", n = x,testing=FALSE)
      expect_no_warning(ssg$get_TransMatrix(old_eqFreqs = old_freqs, new_eqFreqs = new_freqs, testing=TRUE))
=======
      #expect_equal(sum(old_freqs),1)
      #creating a new ssg object to test get_SrMatrix
      ssg <- singleStructureGenerator$new(globalState = "U", n = x,testing=FALSE)
      out <- ssg$get_SrMatrix(old_eqFreqs = old_freqs, new_eqFreqs = new_freqs, testing=TRUE)
>>>>>>> 2da95e729eefa07d7e85da041f38f115f91fe39c
    }
    else{
      test_num <- test_num + 1
    }
  }
})

#testing RE_within
test_that("RE_within test",{
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

test_that("RE_within frequency test",{
  #test whether the eqfrequencies calculated in RE_within make sense by recalculating them afterwards for the FreqsX part
  #case 1
  flag <- TRUE
  Y_seq <- c(rep(1, 1500))
  yEqFreqs <- c(0.9,0.09,0.01)
  while(TRUE){
    ssg <- singleStructureGenerator$new(globalState="U",n=3000)
    out <- ssg$RE_within(Y_seq= Y_seq, Y_eqFreqs = yEqFreqs,testing = TRUE)
    if (out$id == "FreqsY"){
      pi_u <- length(which(ssg$get_seq()[1:1500]==1))/1500
      pi_p <- length(which(ssg$get_seq()[1:1500]==2))/1500
      pi_m <- length(which(ssg$get_seq()[1:1500]==3))/1500
      
      expect_equal(c(pi_u,pi_p,pi_m),yEqFreqs,tolerance = 0.05)
      break
    }
  }
  
  #case 2
  Y_seq_2 <- c(rep(2,500))
  yEqFreqs_2 <- c(0.05,0.9,0.05)
  while(TRUE){
    ssg_2 <- singleStructureGenerator$new(globalState="U",n=2000)
    out_2 <- ssg_2$RE_within(Y_seq= Y_seq_2, Y_eqFreqs = yEqFreqs_2,testing = TRUE)
    
    if(out_2$id =="FreqsY"){
      pi_u_2 <- length(which(ssg_2$get_seq()[1:1500]==1))/1500
      pi_p_2 <- length(which(ssg_2$get_seq()[1:1500]==2))/1500
      pi_m_2 <- length(which(ssg_2$get_seq()[1:1500]==3))/1500
      expect_equal(c(pi_u_2,pi_p_2,pi_m_2),yEqFreqs_2, tolerance=0.05)
      break
    }
  }
  #case 3
  Y_seq_3 <- c(rep(3,100))
  yEqFreqs_3 <- c(0,0.01,0.99)
  while(TRUE){
    ssg_3 <- singleStructureGenerator$new(globalState="U",n=1000)
    out_3 <- ssg_3$RE_within(Y_seq= Y_seq_3, Y_eqFreqs = yEqFreqs_3,testing = TRUE)
    
    if(out_3$id =="FreqsY"){
      pi_u_3 <- length(which(ssg_3$get_seq()[1:900]==1))/900
      pi_p_3 <- length(which(ssg_3$get_seq()[1:900]==2))/900
      pi_m_3 <- length(which(ssg_3$get_seq()[1:900]==3))/900
      expect_equal(c(pi_u_3,pi_p_3,pi_m_3),yEqFreqs_3, tolerance=0.05)
      break
    }
  }
})
  
