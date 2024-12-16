#!/usr/bin/env Rscript
library(optparse)
library(parallel)

#source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\multiRegion_SIM.R")
#source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\nikolas_extenstion.R")
source("/Users/nikolas/Library/CloudStorage/OneDrive-Persönlich/Desktop/Bachelor/MES_GitClone30_10/R/multiRegion_SIM.R")
source("/Users/nikolas/Library/CloudStorage/OneDrive-Persönlich/Desktop/Bachelor/MES_GitClone30_10/R/nikolas_extenstion.R")

option_list <- list(
  make_option("--alpha_pI", type = "numeric", default = 0.1,
              help = "Value of alpha_pI"),
  make_option("--beta_pI", type = "numeric", default = 1,
              help = "Value of beta_pI"),
  make_option("--alpha_mI", type = "numeric", default = 0.1,
              help = "Value of alpha_mI"),  
  make_option("--beta_mI", type = "numeric", default = 0.5,
              help = "Value of beta_mI"),
  make_option("--alpha_pNI", type = "numeric", default = 0.1,
              help = "Value of alpha_pNI"),
  make_option("--beta_pNI", type = "numeric", default = 1,
              help = "Value of beta_pNI"),
  make_option("--alpha_mNI", type = "numeric", default = 0.5,
              help = "Value of alpha_mNI"),
  make_option("--beta_mNI", type = "numeric", default = 0.1,
              help = "Value of beta_mNI"),
  make_option("--alpha_Ri", type = "numeric", default = 0.1,
              help = "Value of alpha_Ri"),
  make_option("--iota", type = "numeric", default = 0.3,
              help = "Value of iota"),
  make_option("--mu", type = "numeric", default = 0,
              help = "Value of mu"),
  make_option("--population", type = "integer", default = NULL,
              help = "Population size"),
  make_option("--samples", type = "integer", default = NULL,
              help = "Sample size"),
  make_option("--nodes", type = "character", default = NULL,
              help = "Nodes directory"),
  make_option("--edges", type = "character", default = NULL,
              help = "Edges directory"),
  make_option("--seqlen", type = "integer", default = NULL,
              help = "Sequence length"),
  make_option("--outputdir", type = "character", default = NULL,
              help = "Output directory"),
  make_option("--nsim", type = "integer", default = NULL,
              help = "Number of simulations"),
  make_option("--ncores", type = "integer", default = NULL,
              help = "Number of cores")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


calculate_mu = function(n_pop, samples_n){
  n_pop <- 0
  ##to be done
}

params <- data.frame(
  alpha_pI = opt[["alpha_pI"]],
  beta_pI = opt[["beta_pI"]],
  alpha_mI = opt[["alpha_mI"]],
  beta_mI = opt[["beta_mI"]],
  alpha_pNI = opt[["alpha_pNI"]],
  beta_pNI = opt[["beta_pNI"]],
  alpha_mNI = opt[["alpha_mNI"]],
  beta_mNI = opt[["beta_mNI"]],
  alpha_Ri =opt[["alpha_Ri"]],
  iota = opt[["iota"]],
  mu =opt[["mu"]]
  
)

i_size <- opt[["seqlen"]] / 5

infoStr <- data.frame(
  n = c(i_size,i_size,i_size,i_size,i_size),
  globalState = c("M","U","M","U","M")  
)


nodes_list <- list.files(path = opt[["nodes"]],full.names = TRUE)
edges_list <- list.files(path = opt[["edges"]],full.names = TRUE)

if(length(nodes_list) != length(edges_list)){
  stop("Nodes and edges files are not matching!")
}

if(!dir.exists(opt[["outputdir"]])){
  dir.create(opt[["outputdir"]],recursive = TRUE)

  }else if(dir.exists(opt[["outputdir"]])){
  unlink(opt[["outputdir"]],recursive = TRUE)
  dir.create(opt[["outputdir"]],recursive = TRUE)
}

pad_n <- nchar(as.character(opt[["nsim"]])) + 1
'
for (i in seq_along(nodes_list)) {
  nodes <- read.csv(nodes_list[[i]])
  edges <- read.csv(edges_list[[i]])
  
  padded_index <- formatC(i, width = pad_n, format = "d", flag = "0")
  file_path <- file.path(opt[["outputdir"]],paste0("nikolas_dataSIM_",padded_index,".RData"))
  
  recombination_handler <- recombMultiRegionSimulator$new(
    nodes = nodes,
    edges = edges,
    infoStr = infoStr,
    params = params,
    dt = 0.01,
    testing = FALSE
  )
  
  tip_list <- recombination_handler$tip_list
  save(tip_list,file = file_path)
}
'


process_simulation <- function(index) {
  if (index > length(nodes_list)) {
    message(sprintf("Index %d exceeds the number of tasks (%d). Skipping.", index, length(nodes_list)))
    return(sprintf("Index %d exceeds the number of tasks (%d). Skipping.", index, length(nodes_list)))
  }
  
  tryCatch({
    nodes <- read.csv(nodes_list[[index]])
    edges <- read.csv(edges_list[[index]])
    
    padded_index <- formatC(index, width = pad_n, format = "d", flag = "0")
    file_path <- file.path(opt[["outputdir"]], paste0("nikolas_dataSIM_", padded_index, ".RData"))
    
    recombination_handler <- recombMultiRegionSimulator$new(
      nodes = nodes,
      edges = edges,
      infoStr = infoStr,
      params = params,
      dt = 0.01,
      testing = FALSE
    )
    
    tip_list <- recombination_handler$tip_list
    save(tip_list, file = file_path)
  }, error = function(e) {
    message(sprintf("Error at index %d: %s", index, e$message))
  })
}

n_cores <- opt[["ncores"]]
mclapply(seq_along(nodes_list),process_simulation,mc.cores = n_cores)



'
n_cores <- opt[["ncores"]]  # Use available cores but not more than tasks
cl <- makeCluster(n_cores)

# Export necessary objects and libraries to the cluster
clusterExport(cl, varlist = c("nodes_list", "edges_list", "pad_n", "opt", "params", "infoStr", "recombMultiRegionSimulator", "process_simulation"))
clusterEvalQ(cl, library(parallel))  # Load parallel package on all nodes

clusterEvalQ(cl, {
  source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\multiRegion_SIM.R")
  source("C:\\Users\\ricca\\OneDrive\\Desktop\\Bachelor\\MES_GitClone30_10\\R\\nikolas_extenstion.R")
})

# Run simulations in parallel
parLapply(cl, seq_along(nodes_list), process_simulation)

# Stop the cluster
stopCluster(cl)


'


























































