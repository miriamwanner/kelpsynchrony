# creating the spline correlogram figure

setwd("/Users/miriam/Desktop/revised_kelp_code/Code") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"

load(file=paste0(resloc, "new_data_processing.RData"))

library(ncdf4)
library(R.matlab)
library(geosphere)
library(plyr)
library(wsyn)
library(ncf)
library(hexbin)
library(ggplot2)
library(gtools)
library(ecodist)
library(PBSmapping)
library(rjson)
library(raptr)
library(ggsn)




# spline correlogram
source("altered_spline_plot_func.R")

pdf(file=paste0(resloc, "new_spline_correlogram.pdf"), width=7, height=6)
spline_corr = Sncf(newLonLatROMS[,1], newLonLatROMS[,2], kelpDataROMSSites, latlon = TRUE) # gets the spline correlogram data
synch_vec <- sm[upper.tri(sm, diag=FALSE)] # this is my synchrony matrix as a vector (for plotting the points)
dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)] # this is my distance matrix as a vector (for plotting the points)
smoothingSpline = smooth.spline(spline_corr$real$predicted$x, spline_corr$real$predicted$y, spar=0.35) # this is the data for the line which I plot after plotting the points
my_spline_plot(spline_corr, ylim=c(-0.25, 1), cex.lab=1.5, cex.axis=1.5, cex.sub=1.5) # this plots the spline correlogram (without points)
points(dist_vec, synch_vec, pch=20, cex=0.2, col=rgb(red = 0, green = 0, blue = 0, alpha=0.1)) # this plots the points on top
lines(smoothingSpline) # this plots the line on top
dev.off()












