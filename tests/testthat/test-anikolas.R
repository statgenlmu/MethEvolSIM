library(devtools)
#testing get_SrMatrix
test_that("get_SrMatrix test",{
test_num <- 1000
  for (x in 1:test_num){
    #creating random vectors with length 3 summing to 1
    random_vector1 <- runif(3) / sum(runif(3))
    random_vector1 <- random_vector1 /sum(random_vector1)
    random_vector2 <- runif(3) / sum(runif(3))
    random_vector2 <- random_vector2 /sum(random_vector2)
    
    old_freqs <- random_vector1
    new_freqs <- random_vector2
    
    if(all.equal(old_freqs, new_freqs)!=TRUE&& sum(old_freqs)==1){
      expect_equal(sum(old_freqs),1)
      #creating a new ssg object to test get_SrMatrix
      ssg <- singleStructureGenerator$new(globalState = "U", n = x,eqFreqs = old_freqs,testing=FALSE)
      out <- ssg$get_SrMatrix(old_eqFreqs = old_freqs, new_eqFreqs = new_freqs, testing=TRUE)
    }
    else{
      test_num <- test_num + 1
    }
  }
})

#testing RE_within
test_that("RE_within test",{
  test_num <- 1000
  for(i in 5:test_num){
    vec <- sample(c(1, 2, 3), i-3, replace = TRUE)
    #creating random vectors summing to 1
    random_vector3 <- runif(3) / sum(runif(3))
    random_vector3 <- random_vector3 /sum(random_vector3)
    random_vector4 <- runif(3) / sum(runif(3))
    random_vector4 <- random_vector4 /sum(random_vector4)
    if(all.equal(random_vector3,random_vector4)!= TRUE && sum(random_vector4)==1){
      ssg <- singleStructureGenerator$new(globalState = "U", n = i,eqFreqs = random_vector4,testing=FALSE)
      newSeq <- ssg$RE_within(Y_seq = vec,Y_eqFreqs = random_vector3,testing=TRUE)
      expect_equal(sum(newSeq$new_obsFreqs), 1)
      
    }
    else{
      test_num <- test_num+1
    }
    
    
  }
 
  
  
  
})
  