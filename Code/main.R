# This script runs the other scripts. Gets the data from "Data" and puts in "Results"
# set working directory of R to this directory

rm(list=ls())

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed

# get the ROMS cells corresponding to the kelp data
source("get_ROMS_from_kelp.R")
source("create_synchrony_matrices.R")
source("intermediate_figures.R") # includes possible figures for the appendix
source("final_figures.R")
source("new_data_processing.R")
