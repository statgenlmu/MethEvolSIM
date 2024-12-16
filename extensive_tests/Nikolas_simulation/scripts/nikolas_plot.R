#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
#library(dplyr)

#library(mgcv)

option_list <- list(
  make_option("--sumstats", type = "character", default = NULL,
              help = "Sumstats file"),
  make_option("--outdir",type = "character", default = NULL,
              help = "Output directory")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sumstats <- read.csv(opt[["sumstats"]])

if(!dir.exists(opt[["outdir"]])){
  dir.create(opt[["outdir"]],recursive = TRUE)
  
}else if(dir.exists(opt[["outdir"]])){
  unlink(opt[["outdir"]],recursive = TRUE)
  dir.create(opt[["outdir"]],recursive = TRUE)
}

y_columns <- setdiff(names(sumstats), c("rec_rate","sim_num","n_rec","tajima_i1","tajima_i2"))

for (y_col in y_columns) {
  name <- ""
  if(y_col == "tajima_island"){
    name <- "Tajimas Pi on Islands"
  }else if (y_col=="tajima_nonisland"){
    name <- "Tajimas Pi on Non-Islands"
  }else if (y_col =="watterson_island_pair"){
    name <- "Tajimas Pi on Pairwise Islands"
  }else if(y_col == "wattersons_island"){
    name <- "Wattersons Theta on Islands"
  }else if(y_col == "wattersons_nonisland"){
    name <- "Wattersons Theta on Non-Islands"
  }
  
  
  plot <- ggplot(sumstats, aes(x = n_rec, y = .data[[y_col]])) +
    # Add jitter only if rec_rate is not constant
    { if (var(sumstats$n_rec, na.rm = TRUE) == 0) {
      geom_point(color = "blue", size = 2, alpha = 0.5)
    } else {
      geom_jitter(width = 0.001, height = 0, color = "blue", size = 2, alpha = 0.5)
    }
    } +
    # Fix y-axis scale from 0 to 1
    #scale_y_continuous(limits = c(0, 1)) +
    # Title and axis labels
    labs(
      title = paste("Recombination Events vs", name),
      x = "Recombination Events",
      y = name
    ) +
    # Use minimal theme with white background
    theme_minimal(base_size = 14) +
    # Customize the theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), # Center title
      axis.title = element_text(face = "bold"), # Bold axis titles
      axis.text = element_text(size = 12), # Larger axis text
      panel.grid.major = element_line(color = "gray", linewidth = 0.5), # Subtle grid
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1) # Add border
    )
  
  

  # Save the plot
  output_file <- file.path(opt[["outdir"]], paste0("scatter_rec_rate_vs_", name, ".png"))
  ggsave(output_file, plot = plot, width = 6, height = 4, dpi = 300)
}

















