# Main script for running experiments to recreate figures
# This script runs the other scripts. Gets the data from "Data" and puts in "Results"

library("here")

here::i_am("kelpsynchrony/Code/main.R")
setwd(here("kelpsynchrony/Code"))
rm(list=ls())
# location for storing the results
resloc <- "../Results/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"

source("init_data_processing.R")
source("spline_figure.R")
source("transparent_plots_figure.R")
source("matrices_figure.R")
source("clustering_figures.R")
source("avg_ts_figures.R")
source("mrm_results.R")
