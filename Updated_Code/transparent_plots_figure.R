# creating the transparent plots figures

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





# plots corrected for other variables

synch_vec <- sm[upper.tri(sm, diag=FALSE)] # this is my synchrony matrix as a vector (for plotting the points)
dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)] # this is my distance matrix as a vector (for plotting the points)


pdf(file=paste0(resloc, "transparent_plots_figure.pdf"), width=7, height=4.2)
par(mfrow=c(2, 3), mai = c(0.4, 0.4, 0.4, 0.1), oma=c(1.5, 1.5, 0, 0))
# layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))

transport_vec <- symProbMat[upper.tri(symProbMat, diag=FALSE)]
log_transport_vec <- logSymProbMat[upper.tri(logSymProbMat, diag=FALSE)]
transport_vec_smooth_spline <- smooth.spline(log_transport_vec, synch_vec, spar=0.8)
# plot(transport_vec, synch_vec, main="Transport", xlab="Transport", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.2))
plot(log_transport_vec, synch_vec, main="Dispersal (log scale)", xlab="log(Transport)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))
lines(transport_vec_smooth_spline, col="#0099FF", lwd=2)

connectivity_vec <- symProbMatF[upper.tri(symProbMatF, diag=FALSE)]
log_connectivity_vec <- logSym[upper.tri(logSym, diag=FALSE)]
connectivity_vec_smooth_spline <- smooth.spline(log_connectivity_vec, synch_vec, spar=0.8)
# plot(connectivity_vec, synch_vec, main="Connectivity", xlab="Connectivity", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.2))
plot(log_connectivity_vec, synch_vec, main="Connectivity (log scale)", xlab="log(Connectivity)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))
lines(connectivity_vec_smooth_spline, col="#0099FF", lwd=2)

dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)]
dist_vec_smooth_spline <- smooth.spline(dist_vec, synch_vec, spar=0.8)
plot(dist_vec, synch_vec, main="Distance (km)", xlab="Distance (km)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))
lines(dist_vec_smooth_spline, col="#0099FF", lwd=2)

waves_vec <- smWaves[upper.tri(smWaves, diag=FALSE)]
waves_vec_smooth_spline <- smooth.spline(waves_vec, synch_vec, spar=0.8)
plot(waves_vec, synch_vec, main="Synchrony in wave height", xlab="Waves", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))
lines(waves_vec_smooth_spline, col="#0099FF", lwd=2)

no3_vec <- smNO3[upper.tri(smNO3, diag=FALSE)]
no3_vec_smooth_spline <- smooth.spline(no3_vec, synch_vec, spar=0.8)
plot(no3_vec, synch_vec, main="Synchrony in nitrate", xlab="Nitrate", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))
lines(no3_vec_smooth_spline, col="#0099FF", lwd=2)

mtext(outer=TRUE, "Predictor", side=1)
mtext(outer=TRUE, "Synchrony", side=2)
dev.off()
