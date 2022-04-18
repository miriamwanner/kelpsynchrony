# script for intermediate figures

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/ROMSKelpData/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"
load(file=paste0(resloc, "create_synchrony.RData"))


library(spatstat)
library(geosphere)
library(wsyn)
library(ecodist)
library(R.matlab)
library(plot3D)
library(RColorBrewer)
library(ggplot2)
library(PBSmapping)
library(plyr)
library(tidyverse)
library(scatterplot3d)
library(rgl)
library(gtools)
library(ncf)
library(ggplot2)
# library(mms)




# THREE BLOCKS FROM SYNCHRONY PLOT - PLOT ON MAP - uses much of Max's code for geography
# first block of sites: from 1 to 15
# second block of sites: from 16 to 20
# last block of sites: from 21 to 43

# part of Max's code for the mapping
# Build map
# Import coastline data from NOAA
limits <- list(x = c(-123 + 360, -114 + 360), y = c(30, 37))
polys <- importGSHHS(gshhsDB="/Users/miriam/Desktop/Data and code for Miriam/gshhg-bin-2.3.7/gshhs_h.b", xlim=limits$x, ylim=limits$y, maxLevel=1, n=0)
polys <- fortify(polys)

border <- importGSHHS(gshhsDB="/Users/miriam/Desktop/Data and code for Miriam/gshhg-bin-2.3.7/wdb_borders_h.b", xlim=limits$x, ylim=limits$y, maxLevel=1, n=0)
border <- fortify(border)

# Plot map
map <- ggplot() + 
  coord_map(xlim=c(-120.8,-117), ylim=c(32.5,34.6), projection="lambert", lat0=32.5, lat1=34.6) +
  geom_polygon(data=polys, aes(x=X, y=Y, group=PID), fill="lightgrey") +
  geom_path(data=polys, aes(x=X, y=Y, group=PID), color="darkgrey", size=0.1) +
  geom_path(data=border, aes(x=X, y=Y, group=PID), color="darkgrey", size=0.1) +
  scale_x_continuous(breaks=c(-117, -118, -119, -120), labels=c("117", "118", "119", "120")) +
  scale_y_continuous(breaks=seq(32.5,34.5, by=0.5))
map

# this is the part I change for the points of what to plot
roms.block <- data.frame(cell=1:43)
roms.block$region <- rep(NA, 43)
roms.block$region[1:15] <- c(rep("block 1", 15)) 
roms.block$region[16:20] <- c(rep("block 2", 5)) 
roms.block$region[21:43] <- c(rep("block 3", 23))
roms.block$i.roms.lat <- newLonLatROMS[,2]
roms.block$i.roms.long <- newLonLatROMS[,1]

# plot the blocks as different colors
roms.map <- map + 
  geom_point(data=roms.block, aes(x=i.roms.long, y=i.roms.lat, fill=region), size=5, shape=21, color="black", alpha=0.5) +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.position=c(0,0), legend.justification=c(0,0), legend.text=element_text(size=18)) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1)) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2)) 
roms.map

# COMPARING SYNCHRONY AT SAME DISTANCE FOR MORE CONNECTED VS LESS CONNECTED SITES variables:
# moreConnected = a vector of all of the values of synchrony for the more connected sites less than 17 km away
# lessConnected = a vector of all of the values of synchrony for the less connected sites less than 17 km away
# moreConnected2 = a vector of all of the values of synchrony for the more connected sites less than 12 km away
# lessConnected2 = a vector of all of the values of synchrony for the less connected sites less than 12 km away


# Want to see if sites that are connected by dispersal are more synchronous than other sites less connected for the same distances.
# plot distance matrix against the probability matrix
x <- as.vector(newDistMat)
y <- as.vector(log10(probMatFecundity))
plot(x[x < 17], y[x < 17]) # just the distances less than 17 km against their probability (top left of previous plot)
# catagory 1 are sites with log10(probability) > -3
# catagory 2 are sites with log10(probability) <= -3
moreConnected <- vector(,51)
lessConnected <- vector(,112)
countMoreConnected <- 0
countLessConnected <- 0
for(i in 1:nrow(newDistMat)){
  for(j in 1:ncol(newDistMat)){
    if(newDistMat[i, j] <= 17){
      if(log10(probMatFecundity[i, j]) > -3){
        countMoreConnected <- countMoreConnected + 1
        moreConnected[countMoreConnected] <- sm[i, j]
      }
      else{
        countLessConnected <- countLessConnected + 1
        lessConnected[countLessConnected] <- sm[i, j]
      }
    }
  }
}
# plot that evenly spaces the points along the x axis
plot(1:length(moreConnected), moreConnected, xlim=c(0,112))
points(1:length(lessConnected), lessConnected, col = "blue")
# plot that stacks the points on top of one another
plot(c(rep(1, 51)) , moreConnected, xlim=c(0,3))
points(c(rep(2, 112)), lessConnected, col = "blue")
boxplot(c(moreConnected,lessConnected) ~ c(rep(1,51),rep(2,112)))

# See if synchrony shows a bigger difference when the sites are closer
plot(x[x < 12], y[x < 12]) # just the distances less than 12 km against their probability (top left of previous plot)
# catagory 1 are sites with log10(probability) > -3
# catagory 2 are sites with log10(probability) <= -3
moreConnected2 <- vector(,40)
lessConnected2 <- vector(,69)
countMoreConnected2 <- 0
countLessConnected2 <- 0
for(i in 1:nrow(newDistMat)){
  for(j in 1:ncol(newDistMat)){
    if(newDistMat[i, j] <= 12){
      if(log10(probMatFecundity[i, j]) > -3){
        countMoreConnected2 <- countMoreConnected2 + 1
        moreConnected2[countMoreConnected2] <- sm[i, j]
      }
      else{
        countLessConnected2 <- countLessConnected2 + 1
        lessConnected2[countLessConnected2] <- sm[i, j]
      }
    }
  }
}
# plot that evenly spaces the points along the x axis
plot(1:length(moreConnected2), moreConnected2, xlim=c(0,70))
points(1:length(lessConnected2), lessConnected2, col = "blue")
# plot that stacks the points on top of one another
plot(c(rep(1, 40)) , moreConnected2, xlim=c(0,3))
points(c(rep(2, 69)), lessConnected2, col = "blue")
boxplot(c(moreConnected2,lessConnected2) ~ c(rep(1,40),rep(2,69)))




# PLOT 3D SCATTER PLOT
# x-axis: distances
# y-axis: log10(symProbMatF)
# z-axis: synchrony

distPlot <- as.vector(newDistMat)
probPlot <- as.vector(logSym)
synchPlot <- as.vector(sm)
nitratePlot <- as.vector(smNO3)
wavesPlot <- as.vector(smWaves)

plot(x=nitratePlot, y=synchPlot)

scatter3D(xDist, yProb, zSync, xlab = "distance", ylab = "symProbMatF", zlab = "synchrony")
scatterplot3d(x=xDist, y=yProb, z=zSync)
options(rgl.printRglwidget = TRUE)
plot3d(x=xDist, y=yProb, z=zSync)
plot3d(x=xDist[xDist < 17], y=yProb[xDist < 17], z=zSync[xDist < 17])


zWaveSync <- as.vector(smWaves)
plot3d(x=xDist, y=yProb, z=zWaveSync)


zNO3Sync <- as.vector(smNO3)
plot3d(x=xDist, y=yProb, z=zNO3Sync)




times <- c(1:ncol(kelpDataROMSSites))

colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
colnames(newLonLatROMS)[2] <- "lat"
coords <- as.data.frame(newLonLatROMS)

cleanKelp <- cleandat(data.matrix(kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp)

kelpClust <- clust(dat=cleanKelp,times=1:matcol,coords=coords,method="pearson")
get_clusters(kelpClust)
plotmap(kelpClust)
# creating wavelet mean fields for each module and plot
kelpClust <- addwmfs(kelpClust)
plotmag(get_wmfs(kelpClust)[[3]][[1]]) # I do not understand the indexing here but I just copied it from the vignette
plotmag(get_wmfs(kelpClust)[[3]][[2]]) 
plotmag(get_wmfs(kelpClust)[[3]][[3]])





# clustering without the sites 41-43
kelpDataROMSSites40 <- kelpDataROMSSites[-c(1, 43), ]
newLonLatROMS40 <- newLonLatROMS[-c(1, 43), ]
times40 <- c(1:ncol(kelpDataROMSSites40))

colnames(newLonLatROMS40)[1] <- "lon" # rename columns for 'clust' method
colnames(newLonLatROMS40)[2] <- "lat"
coords40 <- as.data.frame(newLonLatROMS40)

cleanKelp40 <- cleandat(data.matrix(kelpDataROMSSites40), times40, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp40)

kelpClust40 <- clust(dat=cleanKelp40,times=1:matcol,coords=coords40,method="pearson")
get_clusters(kelpClust40)
plotmap(kelpClust40)
# creating wavelet mean fields for each module and plot
kelpClust40 <- addwmfs(kelpClust40)
plotmag(get_wmfs(kelpClust40)[[3]][[1]]) # I do not understand the indexing here but I just copied it from the vignette
plotmag(get_wmfs(kelpClust40)[[3]][[2]]) 
plotmag(get_wmfs(kelpClust40)[[3]][[3]])





transport <- as.vector(logSymProbMat)
connectivity <- as.vector(logSym)
plot(transport, connectivity) # all the distances against the probability

transport <- as.vector(logSymProbMat)
synchrony <- as.vector(sm)
plot(transport, synchrony)

transport <- as.vector(symProbMat)
synchrony <- as.vector(sm)
plot(transport, synchrony)


hist(as.vector(logSymProbMat), main="Histogram of the Log of Transport")

hist(as.vector(logSym), main="Histogram of the Log of Connectivity")

hist(as.vector(symProbMat), main="Histogram of Transport (not log transformed)")

hist(as.vector(symProbMatF), main="Histogram of Connectivity (not log transformed)")






# plotting synchrony adjusted for nitrate

synch_vec <- sm[upper.tri(sm, diag=FALSE)]
nitrate_vec <-smNO3[upper.tri(smNO3, diag=FALSE)]
transport_vec <- symProbMat[upper.tri(symProbMat, diag=FALSE)]
connectivity_vec <- symProbMatF[upper.tri(symProbMatF, diag=FALSE)]
waves_vec <- smWaves[upper.tri(smWaves, diag=FALSE)]
dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)]

predictors_df_transport <- data.frame(synch_vec, nitrate_vec, transport_vec, waves_vec, dist_vec)
predictors_df_connectivity <- data.frame(synch_vec, nitrate_vec, connectivity_vec, waves_vec, dist_vec)

linearModTransport <- lm(synch_vec ~ nitrate_vec + transport_vec + waves_vec + dist_vec, data=predictors_df_transport)
linearModConnectivity <- lm(synch_vec ~ nitrate_vec + connectivity_vec + waves_vec + dist_vec, data=predictors_df_connectivity)

modTransportCoefficients <- linearModTransport$coefficients
modConnectivityCoefficients <- linearModConnectivity$coefficients


synch_corrected_for_nitrate_transport <- synch_vec - (as.numeric(modTransportCoefficients['nitrate_vec']) * nitrate_vec) - 
  (as.numeric(modTransportCoefficients['waves_vec']) * waves_vec) - (as.numeric(modTransportCoefficients['dist_vec']) * dist_vec)
plot(transport_vec, synch_corrected_for_nitrate_transport, xlab="Transport", ylab="Synchrony corrected for nitrate, waves, distance")
plot(log10(transport_vec), synch_corrected_for_nitrate_transport, xlab="Transport", ylab="Synchrony corrected for nitrate, waves, distance")


synch_corrected_for_nitrate_connectivity <- synch_vec - (as.numeric(modConnectivityCoefficients['nitrate_vec']) * nitrate_vec) - 
  (as.numeric(modTransportCoefficients['waves_vec']) * waves_vec) - (as.numeric(modTransportCoefficients['dist_vec']) * dist_vec)
plot(connectivity_vec, synch_corrected_for_nitrate_connectivity, xlab="Connectivity", ylab="Synchrony corrected for nitrate, waves, distance")
plot(log10(connectivity_vec), synch_corrected_for_nitrate_connectivity, xlab="Connectivity", ylab="Synchrony corrected for nitrate, waves, distance")



# plotting synchrony adjusted for nitrate LOG transforms

log_transport_vec <- logSymProbMat[upper.tri(logSymProbMat, diag=FALSE)]
log_connectivity_vec <- logSym[upper.tri(logSym, diag=FALSE)]

predictors_df_log_transport <- data.frame(synch_vec, nitrate_vec, log_transport_vec, waves_vec, dist_vec)
predictors_df_log_connectivity <- data.frame(synch_vec, nitrate_vec, log_connectivity_vec, waves_vec, dist_vec)

linearModLogTransport <- lm(synch_vec ~ nitrate_vec + log_transport_vec + waves_vec + dist_vec, data=predictors_df_log_transport)
linearModLogConnectivity <- lm(synch_vec ~ nitrate_vec + log_connectivity_vec + waves_vec + dist_vec, data=predictors_df_log_connectivity)

modLogTransportCoefficients <- linearModLogTransport$coefficients
modLogConnectivityCoefficients <- linearModLogConnectivity$coefficients

synch_corrected_for_nitrate_log_transport <- synch_vec - (as.numeric(modTransportCoefficients['nitrate_vec']) * nitrate_vec) - 
  (as.numeric(modTransportCoefficients['waves_vec']) * waves_vec) - (as.numeric(modTransportCoefficients['dist_vec']) * dist_vec)
plot(log_transport_vec, synch_corrected_for_nitrate_log_transport, xlab="Log(Transport)", ylab="Synchrony corrected for nitrate, waves, distance")

synch_corrected_for_nitrate_log_connectivity <- synch_vec - (as.numeric(modTransportCoefficients['nitrate_vec']) * nitrate_vec) - 
  (as.numeric(modTransportCoefficients['waves_vec']) * waves_vec) - (as.numeric(modTransportCoefficients['dist_vec']) * dist_vec)
plot(log_connectivity_vec, synch_corrected_for_nitrate_log_connectivity, xlab="Log(Connectivity)", ylab="Synchrony corrected for nitrate, waves, distance")


save.image(file = paste0(resloc, "intermediate_figures.RData"))

