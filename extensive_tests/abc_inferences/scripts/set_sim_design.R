library(optparse)

# Define command-line options with unique letters
option_list <- list(
  make_option("--dir", type = "character", default = NULL,
              help = "Full path to save the design file", metavar = "character"),
  make_option("--genome-dist", type = "character", default = NULL,
              help = "Distribution of methylation sites in the genomic region", metavar = "character"),
  make_option("--tree", type = "character", default = "all",
              help = "Tree in Newick format", metavar = "character"),
  make_option("--n-sim", type = "integer", default = NULL,
              help = "Number of simulations", metavar = "integer"),
  make_option("--seed", type = "integer", default = NULL,
              help = "Seed value", metavar = "integer"),
  make_option("--alpha_pI_shape", type="numeric", default=3, help="Shape parameter for alpha_pI"),
  make_option("--beta_pI_shape", type="numeric", default=3, help="Shape parameter for beta_pI"),
  make_option("--alpha_pI_scale", type="numeric", default=as.numeric(1/30), help="Scale parameter for alpha_pI"),
  make_option("--beta_pI_scale", type="numeric", default=as.numeric(1/3), help="Scale parameter for beta_pI"),
  make_option("--alpha_pNI_shape", type="numeric", default=3, help="Shape parameter for alpha_pNI"),
  make_option("--beta_pNI_shape", type="numeric", default=3, help="Shape parameter for beta_pNI"),
  make_option("--alpha_pNI_scale", type="numeric", default=as.numeric(1/30), help="Scale parameter for alpha_pNI"),
  make_option("--beta_pNI_scale", type="numeric", default=as.numeric(1/3), help="Scale parameter for beta_pNI"),
  make_option("--alpha_mI_shape", type="numeric", default=3, help="Shape parameter for alpha_mI"),
  make_option("--beta_mI_shape", type="numeric", default=3, help="Shape parameter for beta_mI"),
  make_option("--alpha_mI_scale", type="numeric", default=as.numeric(1/30), help="Scale parameter for alpha_mI"),
  make_option("--beta_mI_scale", type="numeric", default=as.numeric(1/6), help="Scale parameter for beta_mI"),
  make_option("--alpha_mNI_shape", type="numeric", default=3, help="Shape parameter for alpha_mNI"),
  make_option("--beta_mNI_shape", type="numeric", default=3, help="Shape parameter for beta_mNI"),
  make_option("--alpha_mNI_scale", type="numeric", default=as.numeric(1/6), help="Scale parameter for alpha_mNI"),
  make_option("--beta_mNI_scale", type="numeric", default=as.numeric(1/30), help="Scale parameter for beta_mNI"),
  make_option("--alpha_Ri_rate", type="numeric", default=as.numeric(1/2), help="Rate parameter for alpha_Ri"),
  make_option("--mu_rate", type="numeric", default=50, help="Rate parameter for mu"),
  make_option("--iota_min", type="numeric", default=0, help="Minimum value for iota"),
  make_option("--iota_max", type="numeric", default=1, help="Maximum value for iota")
)



# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set the seed
set.seed(opt[["seed"]])

# Generate sampled parameters
sampled_params <- data.frame(
  alpha_pI = rgamma(opt[["n-sim"]], shape = opt[["alpha_pI_shape"]], scale = opt[["alpha_pI_scale"]]), 
  beta_pI = rgamma(opt[["n-sim"]], shape = opt[["beta_pI_shape"]], scale = opt[["beta_pI_scale"]]),
  alpha_pNI = rgamma(opt[["n-sim"]], shape = opt[["alpha_pNI_shape"]], scale = opt[["alpha_pNI_scale"]]),
  beta_pNI = rgamma(opt[["n-sim"]], shape = opt[["beta_pNI_shape"]], scale = opt[["beta_pNI_scale"]]),
  alpha_mI = rgamma(opt[["n-sim"]], shape = opt[["alpha_mI_shape"]], scale = opt[["alpha_mI_scale"]]),
  beta_mI = rgamma(opt[["n-sim"]], shape = opt[["beta_mI_shape"]], scale = opt[["beta_mI_scale"]]),
  alpha_mNI = rgamma(opt[["n-sim"]], shape = opt[["alpha_mNI_shape"]], scale = opt[["alpha_mNI_scale"]]),
  beta_mNI = rgamma(opt[["n-sim"]], shape = opt[["beta_mNI_shape"]], scale = opt[["beta_mNI_scale"]]),
  alpha_Ri = rexp(opt[["n-sim"]], rate = opt[["alpha_Ri_rate"]]),
  mu = rexp(opt[["n-sim"]], rate = opt[["mu_rate"]]),
  iota = runif(opt[["n-sim"]], min = opt[["iota_min"]], max = opt[["iota_max"]])
)

# Get the tree
tree <- opt[["tree"]]

# Load the dataframe containing the genomic distribution of methylation sites
load(file = opt[["genome-dist"]])

# Save the sampled parameters, the tree and the distribution of methylation sites in the design file
save(spatial_str, tree, sampled_params, file = opt[["dir"]])


