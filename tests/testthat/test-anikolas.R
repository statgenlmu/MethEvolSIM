get_private <- function(x){
  x[['.__enclos_env__']]$private
}

#testing get_SrMatrix properties
test_that("get_TransMatrix property test",{
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
      #creating a new ssg object to test get_SrMatrix
      ssg <- singleStructureGenerator$new(globalState = "U", n = x,testing=FALSE)
      expect_no_warning(ssg$get_TransMatrix(old_eqFreqs = old_freqs, new_eqFreqs = new_freqs, testing=TRUE))
    }
    else{
      test_num <- test_num + 1
    }
  }
})

#testing RE_within properties with random vectors
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
test_that("RE_within frequency test",{

  #different testing cases
  testing_parameters <- list(list(Y_seq = c(rep(1, 1500)), yEqFreqs = c(0.9,0.09,0.01), globalState = "U",n=3000),
                             list(Y_seq = c(rep(2,500)), yEqFreqs = c(0.05,0.9,0.05), globalState = "U",n=2000),
                             list(Y_seq = c(rep(3,100)), yEqFreqs = c(0,0.01,0.99), globalState = "U",n=1000),
                             list(Y_seq = c(rep(2,200)), yEqFreqs = c(0,0.99,0.01), globalState = "M",n=400)
                             )

  for(i in 1:length(testing_parameters)){
    Y_seq_test <- testing_parameters[[i]]$Y_seq
    yEqFreqs_test <- testing_parameters[[i]]$yEqFreqs
    globalState_test <- testing_parameters[[i]]$globalState
    n_test  <- testing_parameters[[i]]$n
    ssg <- singleStructureGenerator$new(globalState=globalState_test,n=n_test)
    xEqFreqs_test <- ssg$get_eqFreqs()
    out <- ssg$RE_within(Y_seq= Y_seq_test, Y_eqFreqs = yEqFreqs_test,testing = TRUE)
    
    #when the xfreqs were chosen we calculate the expected frequencies for the yseq part
    if (out$id =="FreqsX"){
      xSeq <- out$chosenSeq
      len_x <- length(xSeq)
      len_y <- length(ssg$get_seq())-len_x
      
      pi_u_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==1))/len_y
      pi_p_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==2))/len_y
      pi_m_x <- length(which(ssg$get_seq()[len_x+1:length(ssg$get_seq())]==3))/len_y
      expect_equal(c(pi_u_x,pi_p_x,pi_m_x),xEqFreqs_test,tolerance=0.2)

      #X_seq_in <- out$chosenSeq
      #X_seq_out <- out$newSeq[1:length(X_seq_in)]
      #expect_equal(X_seq_in,X_seq_out)
    }
    #when the yfreqs were chosen we calculate the expected frequencies for the xseq part
    if (out$id == "FreqsY"){
      len <- length(Y_seq_test)
      pi_u <- length(which(ssg$get_seq()[1:len]==1))/len
      pi_p <- length(which(ssg$get_seq()[1:len]==2))/len
      pi_m <- length(which(ssg$get_seq()[1:len]==3))/len
      expect_equal(c(pi_u,pi_p,pi_m),yEqFreqs_test,tolerance = 0.2)
    }
  }
})

#here i test wether the X and Y sequence that is not part of the recombining ssg stays the same after recombination as it should
test_that("RE_complete unchanged X/Y sequence test",{
  #testing cases
  testing_parameters <- list(list(infoStr_1 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),position = 30),
                             list(infoStr_1 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),infoStr_2 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),position = 3),
                             list(infoStr_1 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),position = 3),
                             list(infoStr_1 = data.frame(n = c(1000,10),globalState = c("U","U")),infoStr_2 = data.frame(n = c(1000,10),globalState = c("M","M")),position = 1004)
                            )
  
  for(i in 1:length(testing_parameters)){
    #testing for the Xseq
    infoStr_1_testing <- testing_parameters[[i]]$infoStr_1
    infoStr_2_testing <- testing_parameters[[i]]$infoStr_2
    position_testing <- testing_parameters[[i]]$position
    
    csg_1 <- combiStructureGenerator$new(infoStr = infoStr_1_testing)
    csg_1_copy <- csg_1$copy()
    csg_2 <- combiStructureGenerator$new(infoStr = infoStr_2_testing)
    csg_2_copy <- csg_2$copy()

    out <- csg_1$RE_complete(combi_2 = csg_2,position = position_testing,testing=TRUE)
    ssg_position <- out$ssg_index
    ssg_within_position <- out$ssg_position_within
    
    seq_before_x <- c()
    seq_after_x <- c()

    eqFreqs_before_x <- c()
    eqFreqs_after_x <- c()

    for(i in 1:ssg_position){
      if(i == ssg_position){
        #i think i can skip this since this is already tested in the RE_within test.
        next
        #these two lines down below are used for appending the sequence in the recombining ssg but id need to include additional
        #information about which sequence got chosen during the RE_within process and that was already tested above so i skipped it
        #update: it probably wouldnt make sense since the neighborstates get updated and thus the sequences differ from the earlier ones
        #seq_before_x <- append(seq_before_x, csg_1_copy$get_singleStr(i)$get_seq()[1:ssg_within_position-1])
        #seq_after_x <- append(seq_after_x,csg_1$get_singleStr(i)$get_seq()[1:ssg_within_position-1])
      }
      else{
        seq_before_x <- append(seq_before_x, csg_1_copy$get_singleStr(i)$get_seq())
        seq_after_x <- append(seq_after_x, csg_1$get_singleStr(i)$get_seq())
        
        eqFreqs_before_x <- append(eqFreqs_before_x, csg_1_copy$get_singleStr(i)$get_eqFreqs())
        eqFreqs_after_x <- append(eqFreqs_after_x, csg_1$get_singleStr(i)$get_eqFreqs())
      }
    }
    
    expect_equal(seq_before_x,seq_after_x,info="unaffected X sequence is not equal after recombination")
    expect_equal(eqFreqs_before_x,eqFreqs_after_x,info="unaffected X eqFreqs are not equal after recombination")
  #testing for the Yseq
  seq_before_y <- c()
  seq_after_y <- c()
  
  eqFreqs_before_y <- c()
  eqFreqs_after_y <- c()
  
  for(i in ssg_position:csg_1$get_singleStr_number()){
    if(i==ssg_position){
      #i think i can skip this since this is already tested in the RE_within test.
      next
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

#testing wether the recombining ssg has the correct equilibrium frequecies
test_that("RE_complete eqFreqs test",{
  testing_parameters <- list(list(infoStr_1 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4,7,8),globalState = c("M","M","U","M","M","U","M")),position = 30),
                             list(infoStr_1 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),infoStr_2 = data.frame(n = c(2,11,2),globalState = c("M","M","U")),position = 3),
                             list(infoStr_1 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),infoStr_2 = data.frame(n = c(13,5,7,3,4),globalState = c("M","M","U","M","M")),position = 3),
                             list(infoStr_1 = data.frame(n = c(1000,10),globalState = c("U","U")),infoStr_2 = data.frame(n = c(1000,10),globalState = c("M","M")),position = 1004)
  )
  
  for(i in 1:length(testing_parameters)){
    infoStr_1_testing <- testing_parameters[[i]]$infoStr_1
    infoStr_2_testing <- testing_parameters[[i]]$infoStr_2
    position_testing <- testing_parameters[[i]]$position
    
    csg_1 <- combiStructureGenerator$new(infoStr = infoStr_1_testing)
    csg_2 <- combiStructureGenerator$new(infoStr = infoStr_2_testing)

    out <- csg_1$RE_complete(combi_2 = csg_2,position = position_testing,testing=TRUE)
    ssg_position <- out$ssg_index
    ssg_eqFreqs <- out$recombEqFreqs
    
    expect_equal(ssg_eqFreqs, csg_1$get_singleStr(ssg_position)$get_eqFreqs(),info="EqFreqs are not correctly applied to recombining ssg")
  }
})


#test whether the number of outputs of recombMultiRegionSimulator is equal to nodes with time 0 in the msprime dataframe
test_that("recombMultiRegionSimulator output number test",{
  #msprime tables and input parameters
  nodes_path <- system.file("extdata","nodes1_test_17.csv",package = "MethEvolSIM")
  nodes <- read.csv(nodes_path)
  edges_path <- system.file("extdata","edges1_test_17.csv",package = "MethEvolSIM")
  edges <- read.csv(edges_path)
  
  expected_num <- sum(nodes$time == 0)
  unique(nodes$flags)
  infoStr <- data.frame(
    n = c(3,3),
    globalState = c("M","U")
  )
  
  params <- data.frame(
    alpha_pI = 0.1,
    beta_pI = 1,
    alpha_mI = 0.1,
    beta_mI = 0.5,
    alpha_pNI = 0.1,
    beta_pNI = 1,
    alpha_mNI = 0.5,
    beta_mNI = 0.1,
    alpha_Ri = 0.1,
    iota = 0.3,
    mu = 0.1
  )
  
  recombination_handler <- recombMultiRegionSimulator$new(
    nodes = nodes,
    edges = edges,
    infoStr = infoStr,
    params = params,
    dt = 0.01,
    testing = TRUE
  )
  
  result_list <- recombination_handler$tip_list
  #number of nodes with time 0 needs to match the number of nodes in the output
  expect_equal(expected_num, length(result_list))
})
  
#test whether the order of recombination events is correct. also tests for correct order of tips and non tips.
test_that("recombMultiRegionSimulator order of events test",{
  
  nodes_path <- system.file("extdata","nodes1_test_18.csv",package = "MethEvolSIM")
  nodes <- read.csv(nodes_path)
  edges_path <- system.file("extdata","edges1_test_18.csv",package = "MethEvolSIM")
  edges <- read.csv(edges_path)
  
  #i assume that the nodes file always gets output in ascending order by in the ID column. otherwise the test would fail
  
  infoStr <- data.frame(
    n = c(3,3),
    globalState = c("M","U")
  )
  
  params <- data.frame(
    alpha_pI = 0.1,
    beta_pI = 1,
    alpha_mI = 0.1,
    beta_mI = 0.5,
    alpha_pNI = 0.1,
    beta_pNI = 1,
    alpha_mNI = 0.5,
    beta_mNI = 0.1,
    alpha_Ri = 0.1,
    iota = 0.3,
    mu = 0.1
  )
  
  recombination_handler <- recombMultiRegionSimulator$new(
    nodes = nodes,
    edges = edges,
    infoStr = infoStr,
    params = params,
    dt = 0.01,
    testing = TRUE
  )
  
  order_list <- recombination_handler$testing_list
  expect_equal(order_list$order,as.list(nodes$flags))

})


#test if number of ssgs in all processed csgs is as expected
test_that("recombMultiRegionSimulator seq len test",{
  nodes_path <- system.file("extdata","nodes1_test_18.csv",package = "MethEvolSIM")
  nodes <- read.csv(nodes_path)
  edges_path <- system.file("extdata","edges1_test_18.csv",package = "MethEvolSIM")
  edges <- read.csv(edges_path)
  
  #i assume that the nodes file always gets output in ascending order by in the ID column. otherwise the test would fail
  
  infoStr <- data.frame(
    n = c(3,3),
    globalState = c("M","U")
  )
  
  params <- data.frame(
    alpha_pI = 0.1,
    beta_pI = 1,
    alpha_mI = 0.1,
    beta_mI = 0.5,
    alpha_pNI = 0.1,
    beta_pNI = 1,
    alpha_mNI = 0.5,
    beta_mNI = 0.1,
    alpha_Ri = 0.1,
    iota = 0.3,
    mu = 0.1
  )
  
  recombination_handler <- recombMultiRegionSimulator$new(
    nodes = nodes,
    edges = edges,
    infoStr = infoStr,
    params = params,
    dt = 0.01,
    testing = TRUE
  )
  
  node_list <- recombination_handler$node_list
  for(node in node_list){
    #test whether the numbers of ssg in combis are as expected
    ssg_num <- node$combi$get_singleStr_number()
    expect_equal(2,ssg_num,info = "Ssg numbers in combi are not as expected!")
    
    #check whether the length of the ssgs are as expected
    for(j in ssg_num){
      ssgs <- node$combi$get_singleStr(j)
      expect_equal(length(ssgs$get_seq()),3)
    }
    
    #test whether the global states in all processed nodes are as expected.
    #they shouldnt be changed during the process
    global_state <- get_private(node$combi)$singleStr_globalState
    expect_equal(global_state, infoStr$globalState, info="Global states are not as expected!")
    
  }
  
 
})

#test wether all single structures in recombining combistructures do have the same id as their combistructure
test_that("single structure and combistructure id test",{
  #sample file
  nodes_path <- system.file("extdata","nodes3.csv",package = "MethEvolSIM")
  nodes <- read.csv(nodes_path)
  edges_path <- system.file("extdata","edges3.csv",package = "MethEvolSIM")
  edges <- read.csv(edges_path)
  
  infoStr <- data.frame(
    n = c(15,15),
    #u is for island m is for unmethylated i need 5 -> mumum i choose n
    globalState = c("M","U")
  )
  
  params <- data.frame(
    alpha_pI = 0.1,
    beta_pI = 1,
    alpha_mI = 0.1,
    beta_mI = 0.5,
    alpha_pNI = 0.1,
    beta_pNI = 1,
    alpha_mNI = 0.5,
    beta_mNI = 0.1,
    alpha_Ri = 0.1,
    iota = 0.3,
    mu = 0.1
  )
  
  recombination_handler <- recombMultiRegionSimulator$new(
    nodes = nodes,
    edges = edges,
    infoStr = infoStr,
    params = params,
    dt = 0.01,
    testing = TRUE
  )
  testing_list <- recombination_handler$testing_list$combis$combi_1
  print(paste0("Length testing list: ", length(testing_list)))
  
  for(index in 1:length(testing_list)){
    combi_id <- testing_list[[index]]$get_id()
    ssg_num <- testing_list[[index]]$get_singleStr_number()
    #print(paste0("combi_id ", combi_id, " ssg_num ", ssg_num))
    for(k in 1:ssg_num){
      ssg_id <- testing_list[[index]]$get_singleStr(k)$get_myCombiID()
      print(ssg_id)
      expect_equal(combi_id,ssg_id)
      
    }
  }
})















  

