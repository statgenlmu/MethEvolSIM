#!/usr/bin/env Rscript
library(optparse)
#source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\multiRegion_SIM.R")
#source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\nikolas_extenstion.R")
source("/Users/nikolas/Library/CloudStorage/OneDrive-Persönlich/Desktop/Bachelor/MES_GitClone30_10/R/multiRegion_SIM.R")
source("/Users/nikolas/Library/CloudStorage/OneDrive-Persönlich/Desktop/Bachelor/MES_GitClone30_10/R/nikolas_extenstion.R")

option_list <- list(
  make_option("--tipdir", type = "character", default = NULL,
              help = "Dir of tips list"),
  make_option("--output", type = "character", default = NULL,
              help = "Dir of output"),
  make_option("--nsim", type = "integer", default = NULL,
              help = "Number of simulations"),
  make_option("--recrates", type = "character", default = NULL,
              help = "Recrates location"),
  make_option("--samples", type = "integer", default = NULL,
              help = "Sample number")
  )

convert_global_methylation = function(seq){
  seq[seq==1] <- 0
  seq[seq==2] <- 0.5
  seq[seq==3] <- 1
  
  return(seq)
}
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
tiplist <- list.files(path = opt[["tipdir"]],full.names = TRUE)

rec_rate_df <- read.csv(opt[["recrates"]])

samples <- opt[["samples"]]
nsim <- opt[["nsim"]]
pad_n <- nchar(as.character(length(tiplist))) + 1

#results dataframe
results <- data.frame(sim_num = integer(), rec_rate = numeric(),n_rec = numeric(), tajima_island = numeric(), tajima_non_island = numeric(), wattersons_island = numeric(), wattersons_nonisland = numeric(), watterson_island_pair = numeric(), tajima_i1 = numeric(), tajima_i2 = numeric(), stringsAsFactors = FALSE)

tips_list <- list.files(path = opt[["tipdir"]],full.names = TRUE)


compute_tajima <- function(tips_list,tip_num){
  #list the labels o or u, over 0.5 or under 0.5. for each simulation
  #loop each tip_list
  for(i in seq_along(tips_list)){
    label_list_islands <- list(rep(NA,tip_num),rep(NA,tip_num))
    label_list_nonislands <- list(rep(NA,tip_num),rep(NA,tip_num),rep(NA,tip_num))
    
    test_list_islands <- list(rep(NA,tip_num),rep(NA,tip_num))
    
    load(tips_list[[i]])
    sim_num <- sub(".*_([0-9]+)\\.RData$", "\\1", tips_list[[i]])
    nodes <- read.csv(paste0("./out/nodes/nodes_",sim_num,".csv"))
    n_rec <- sum(nodes$flags == 131072, na.rm = TRUE)/2
    #loop each tip
    for(j in 1:length(tip_list)){
      #loop each ssg
      for(k in 1:tip_list[[j]]$combi$get_singleStr_number()){
        if(k == 2 || k==4){
          
          if(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))<0.5){
            label_list_islands[[k/2]][[j]] <- "u"
            test_list_islands[[k/2]][[j]] <- mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))
          }
          else{
            label_list_islands[[k/2]][[j]] <- "o"
            test_list_islands[[k/2]][[j]] <- mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))
          }
        }
        else if(k==1){
          if(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))<0.5){
            label_list_nonislands[[1]][[j]] <- "u"
          }
          else{
            label_list_nonislands[[1]][[j]] <- "o"
          }
        }
        else if(k==3){
          if(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))<0.5){
            label_list_nonislands[[2]][[j]] <- "u"
          }
          else{
            label_list_nonislands[[2]][[j]] <- "o"
          }
        }
        else if(k==5){
          if(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(k)$get_seq()))<0.5){
            label_list_nonislands[[3]][[j]] <- "u"
          }
          else{
            label_list_nonislands[[3]][[j]] <- "o"
          }
        }
      }
    }
    #tajimas pi
    '
    i1_list <- list()
    i2_list <- list()
    for (i in 1:length(label_list_islands)) {
      i1_list <- append(i1_list,label_list_islands[[i]][[1]])
      i2_list <- append(i2_list,label_list_islands[[i]][[2]])
    }
    
    ni1_list <- list()
    ni2_list <- list()
    ni3_list <- list()
    for (i in 1:length(label_list_nonislands)) {
      ni1_list <- append(ni1_list,label_list_nonislands[[i]][[1]])
      ni2_list <- append(ni2_list,label_list_nonislands[[i]][[2]])
      ni3_list <- append(ni3_list,label_list_nonislands[[i]][[3]])
    }
    '
    
    #tajima_island_z <- sum(unlist(i1_list) == "u") * sum(unlist(i1_list) == "o") + sum(unlist(i2_list) == "u") * sum(unlist(i2_list) == "o")
    tajima_i1_z <- sum(unlist(label_list_islands[[1]]) == "u")* sum(unlist(label_list_islands[[1]]) == "o")
    tajima_i2_z <- sum(unlist(label_list_islands[[2]]) == "u")* sum(unlist(label_list_islands[[2]]) == "o")
    
    tajima_i2_n <- choose(samples,2)
    tajima_i1 <- tajima_i1_z/tajima_i2_n
    tajima_i2 <- tajima_i2_z/tajima_i2_n
    
    tajima_island_z <- sum(unlist(label_list_islands[[1]]) == "u")* sum(unlist(label_list_islands[[1]]) == "o") + sum(unlist(label_list_islands[[2]]) == "u")* sum(unlist(label_list_islands[[2]]) == "o")
    tajima_island_n <- choose(samples,2)*2
    tajima_island <- tajima_island_z/tajima_island_n
    
    tajima_nonisland_z <- sum(unlist(label_list_nonislands[[1]]) == "u")* sum(unlist(label_list_nonislands[[1]]) == "o")+sum(unlist(label_list_nonislands[[2]]) == "u")* sum(unlist(label_list_nonislands[[2]]) == "o")+sum(unlist(label_list_nonislands[[3]]) == "u")* sum(unlist(label_list_nonislands[[3]]) == "o")
    tajima_nonisland_n <- choose(samples,2)*3
    tajima_nonisland <- tajima_nonisland_z/tajima_nonisland_n
    
    
    #wattersons tetha
    wattersons_island_s_2 <- if (length(unique(unlist(label_list_islands[[1]]))) == 1) 0 else 1
    wattersons_island_s_4 <- if (length(unique(unlist(label_list_islands[[2]]))) == 1) 0 else 1
    wattersons_island_z <- wattersons_island_s_2 + wattersons_island_s_4
    
    wattersons_nonisland_s_1 <- if (length(unique(unlist(label_list_nonislands[[1]]))) == 1) 0 else 1
    wattersons_nonisland_s_3 <- if (length(unique(unlist(label_list_nonislands[[2]]))) == 1) 0 else 1
    wattersons_nonisland_s_5 <- if (length(unique(unlist(label_list_nonislands[[3]]))) == 1) 0 else 1
    wattersons_nonisland_z <- wattersons_nonisland_s_1 +wattersons_nonisland_s_3 + wattersons_nonisland_s_5
    
    harmonic_sum <- 0
    for(i in 1:(tip_num-1)){
      harmonic_sum <- harmonic_sum + 1/i
    }
    
    wattersons_island <- wattersons_island_z/harmonic_sum
    wattersons_nonisland <- wattersons_nonisland_z/harmonic_sum
    
    #wattersons pairwise islands
    island_pairlist <- mapply(paste0, label_list_islands[[1]],label_list_islands[[2]])
    watterson_island_pair_z <- 0
    for(e in 1:(length(island_pairlist)-1)){
      for(f in (e+1):length(island_pairlist)){
        if(island_pairlist[[e]] != island_pairlist[[f]]){
          watterson_island_pair_z <- watterson_island_pair_z +1
        }
      }
    }
    
      'sum(unlist(island_pairlist)=="uu")*(length(island_pairlist) - sum(unlist(island_pairlist)=="uu")) + 
      sum(unlist(island_pairlist)=="ou")*(length(island_pairlist) - sum(unlist(island_pairlist)=="ou")) +
      sum(unlist(island_pairlist)=="oo")*(length(island_pairlist) - sum(unlist(island_pairlist)=="oo")) +
      sum(unlist(island_pairlist)=="uo")*(length(island_pairlist) - sum(unlist(island_pairlist)=="uo"))'
    
    watterson_island_pair_n <- choose(tip_num,2)
    watterson_island_pair <- watterson_island_pair_z/watterson_island_pair_n
    
    
    rec_rate <- rec_rate_df$recrate[rec_rate_df$sim_id == as.integer(sim_num)]
    results <- rbind(results, data.frame(sim_num = sim_num, rec_rate = rec_rate,n_rec = n_rec, tajima_island = tajima_island, tajima_nonisland = tajima_nonisland,wattersons_island = wattersons_island, wattersons_nonisland = wattersons_nonisland,watterson_island_pair = watterson_island_pair, tajima_i1=tajima_i1,tajima_i2=tajima_i2))
    
    if(tajima_island != 0 || tajima_island == 0){
      print("------------------------------------------")
      print(paste0("Tajimas pi: ", tajima_island))
      for(j in 1:length(tip_list)){
        print(paste0("Equilibrium Frequencies of Island 2 at tip ", j))
        print(tip_list[[j]]$combi$get_singleStr(2)$get_eqFreqs())
        print(" ")
        print(paste0("Equilibrium Frequencies of Island 4 at tip ", j))
        print(tip_list[[j]]$combi$get_singleStr(4)$get_eqFreqs())
      }
      print("Global Methylation Stats of each tip: ")
      print(test_list_islands)
      #print(paste0("Tajima island is ", tajima_island))
      #print(label_list_islands)
    }
  }
  return(results)
}

tajima_res <- compute_tajima(tips_list =tips_list, tip_num = samples )

write.csv(tajima_res, paste0(opt[["output"]],"/summary_stats.csv"), row.names = FALSE)
cat(paste0("Results have been saved to ",opt[["output"]],".csv"))




















'
for (i in seq_along(tips_list)) {
  load(tips_list[[i]])
  sim_num <- sub(".*_([0-9]+)\\.RData$", "\\1", tips_list[[i]])
  diff <- 0
  iterations <-0
  padded_index <- formatC(i, width = pad_n, format = "d", flag = "0")
  for (i in 1:(length(tip_list)-1)) {
    tip_1_islands <- list(mean(convert_global_methylation(tip_list[[i]]$combi$get_singleStr(2)$get_seq())),mean(convert_global_methylation(tip_list[[i]]$combi$get_singleStr(4)$get_seq())))
    tip_1_nonislands <-list(mean(convert_global_methylation(tip_list[[i]]$combi$get_singleStr(1)$get_seq())),mean(convert_global_methylation(tip_list[[i]]$combi$get_singleStr(3)$get_seq())),mean(convert_global_methylation(tip_list[[i]]$combi$get_singleStr(5)$get_seq())))
    for(j in (i+1):length(tip_list)){
      iterations <- iterations+1
      #list of two elements representing both islands umumu
      tip_2_islands <- list(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(2)$get_seq())),mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(4)$get_seq())))
      tip_2_nonislands <-list(mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(1)$get_seq())),mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(3)$get_seq())),mean(convert_global_methylation(tip_list[[j]]$combi$get_singleStr(5)$get_seq())))
      
      
      
      
      
      for(k in length(tip_1_islands)){
        if(tip_1_islands[k]>0.5 && tip_2_islands[k]<0.5){
          diff <- diff+1
          break
        }
        else if(tip_1_islands[k]<0.5 && tip_2_islands[k]>0.5){
          diff <- diff+1
          break
        }
      }
      
    }
  
  }
  
  
  tajima <- diff/iterations
  rec_rate <- rec_rate_df$recrate[rec_rate_df$sim_id == as.integer(sim_num)]
  results <- rbind(results, data.frame(sim_num = padded_index, rec_rate = rec_rate, tajima = tajima))
  
  
}
'





























































