# script for creating the synchrony matrices

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/ROMSKelpData/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"
load(file=paste0(resloc, "get_ROMS.RData"))

survival_rate <- 0.1

library(spatstat)
library(geosphere)
library(spam)
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



# ==========================================================================================
# Doing stuff with the data

# FECUNDITY AND PROBABILITY MATRIX variables:
# sqrtKelp = square root of all the cells from the summed kelp/ROMS sites
# kelpFecundity = vector of each time series averaged and multiplied by 1463 - each cell is one site
# probMat = each cell is (0.1)^(days in connectivity matrix), diagonal is zero
# probMatFecundity = the probability matrix where each element is multiplied by the fecundity of each donor patch (similar to formula 1 in table of paper)

# square root of every element in the biomass time series for each location
sqrtKelp <- kelpDataROMSSites
for(i in 1:nrow(kelpDataROMSSites)){
  for(j in 1:ncol(kelpDataROMSSites)){ # start at 3 because first column is the site numbers 
    # and take difference from adjacent points
    sqrtKelp[i, j] = sqrt(kelpDataROMSSites[i, j])
  }
}
# average the time series and multiply by 1463 to get a vector of the location's fecundity 
# (I think this is the fecundity - figure 2 of paper)
kelpFecundity <- 1463 * (rowSums(sqrtKelp) / ncol(sqrtKelp))

# create a probability matrix and probability matrix multiplied by fecundity
probMat <- (survival_rate)^oceanAvgKelpSites
diag(probMat) <- 0
# multiply each element of the probability matrix by the fecundity of each donor patch
# (similar to formula 1 in table of paper)
probMatFecundity <- probMat
for(i in 1:nrow(probMatFecundity)){
  probMatFecundity[i,] <- (probMat[i,] * kelpFecundity[i])
}

# DISTANCES FROM ONE CELL TO ANOTHER variables:
# lonLatROMS = matrix (two columns) of all the ROMS sites longitude (column 1) and latitude (column 2)
# newLonLatROMS = matrix (two columns) of the ROMS sites being used for kelp data longitude (column 1) and latitude (column 2)
# newDistMat = the distance matrix of the sites being used for longitude and latitude

# find all distance from cells to each other in a matrix (will be symmetric)
lonLatROMS <- matrix(c(as.vector(unlist(romsCoordinates$lon))[1:62], as.vector(unlist(romsCoordinates$lat))[1:62]), ncol = 2)
newLonLatROMS <- lonLatROMS[1:43,]
count <- 1
for(i in 1:nrow(lonLatROMS)){
  if(is.element(i, newSites)){
    newLonLatROMS[count,] <- lonLatROMS[i,]
    count <- count + 1
  }
}
newDistMat <- distm(newLonLatROMS, newLonLatROMS, fun = distHaversine) * (0.001) # multiply to get km

# plot distance matrix against the probability matrix 
x <- as.vector(newDistMat)
y <- as.vector(log10(probMatFecundity))
plot(x, y) # all the distances against the probability
plot(x[x < 50], y[x < 50]) # just the distances less than 50 km against their probability (top left of previous plot)

# SYNCHRONY MATRIX variables:
# kelpDataROMSSites = the kelp data but without the 'Site' column
# times = the times that will be passed into the cleandat function
# matKelp = the data fram kelpDataROMSSites as a matrix
# cleanKelp = the kelp data cleaned by cleandat() with a cleaning level of two
# matcol = the number of columns in cleanKelp
# sm = the synchrony matrix from using synmat()

# process the data for the synchrony matrix
kelpDataROMSSites <- subset(kelpDataROMSSites, select=-Site) # remove the 'Site' column
times <- c(1:ncol(kelpDataROMSSites)) # this is the times that will be passed into cleandat
matKelp <- data.matrix(kelpDataROMSSites) # creates the matrix from the dataframe
cleanKelp <- cleandat(matKelp, times, clev=2)$cdat # cleans the data - 2 is enough
matcol <- ncol(cleanKelp)
sm <- synmat(cleanKelp, 1:matcol, method="pearson")
diag(sm) <- 1
fields::image.plot(1:ncol(sm),1:ncol(sm),sm,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)

# SYMMETRIC PROBABILITY*FECUNDITY MATRIX variables:
# symProbMatF = a symmetric matrix of probMatFecundity which just takes the max of [i, j] and [j, i] from probMatFecundity
# logSym = the log10 of the symProbMatF matrix - diagonal is 1

# now I make a symmetric matrix of the probMatFecundity matrix by just using the max of [i, j] and [j, i]
symProbMatF <- pmax(probMatFecundity, t(probMatFecundity))
logSym <- log10(symProbMatF)
diag(logSym) <- 1
fields::image.plot(1:ncol(symProbMatF),1:ncol(symProbMatF),logSym,col=heat.colors(20))
fields::image.plot(1:ncol(logSym),1:ncol(logSym),logSym,col=heat.colors(20))




# SYMMETRIC PROBABILITY MATRIX variables:
# symProbMat = a symmetric matrix of probMat matrix by just using the max of [i, j] and [j, i] from probMat
# logSymProbMat = the log10 of the symProbMat matrix - diagonal is 1
symProbMat <- pmax(probMat, t(probMat))
logSymProbMat <- log10(symProbMat)
diag(logSym) <- 1

# WAVES DATA variables:
# waves = the original waves data
# wavesPC = wave data ordered by paddle coordinate
# wavesROMSSites = wave data that matches up to the ROMS sites

# read in the waves data
waves <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Hs_Max_2019.csv")
wavesPC <- waves[match(order.paddle, waves$Site),]
wavesPC <- wavesPC[-whichSitesNA]

# get the waves sites to line up with the ROMS sites (same way as kelp sites)
wavesROMSSites <- wavesPC[0:43,]
counter <- 0
for(i in 1:nrow(wavesPC)){
  loc <- which(vecLocations == i)
  if(length(loc) != 0){
    counter <- counter + 1
    wavesROMSSites[counter,] <- colSums(x = wavesPC[loc,], na.rm = TRUE) / length(loc) 
    # should I average or add for waves? I think average because it is with height
  }
}

# WAVES SYNCHRONY MATRIX variables:
# timesWaves = the times that will be passed into cleandat
# matWaves = the matrix created from the dataframe wavesROMSSites
# cleanWaves = the cleaned wave data using cleandat
# waveCol = the number of columns in the clean wave data
# smWaves = the synchrony matrix created from the clean wave data

# wave data for the synchrony matrix of wave data
wavesROMSSites <- subset(wavesROMSSites, select=-Site) # remove the 'Site' column
timesWaves <- c(1:ncol(wavesROMSSites)) # this is the times that will be passed into cleandat
matWaves <- data.matrix(wavesROMSSites) # creates the matrix from the dataframe
cleanWaves <- cleandat(matWaves, timesWaves, clev=2)$cdat # cleans the data - 2 is enough
waveCol <- ncol(cleanWaves)
smWaves <- synmat(cleanWaves, 1:waveCol, method="pearson")
diag(smWaves) <- 1
fields::image.plot(1:ncol(smWaves),1:ncol(smWaves),smWaves,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)


fields::image.plot(1:ncol(symProbMat),1:ncol(symProbMat),symProbMat,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)
fields::image.plot(1:ncol(logSym),1:ncol(logSym),logSym,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)



# NITRATE DATA variables:
# NO3 = the original nitrate data
# NO3PC = nitrate data ordered by paddle coordinate
# NO3ROMSSites = nitrate data that matches up to the ROMS sites

# read in the waves data
NO3 <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/NO3_Mean_2019.csv")
NO3PC <- NO3[match(order.paddle, NO3$Site),]
NO3PC <- NO3PC[-whichSitesNA]

# get the nitrate sites to line up with the ROMS sites (same way as kelp sites)
NO3ROMSSites <- NO3PC[0:43,]
counter <- 0
for(i in 1:nrow(NO3PC)){
  loc <- which(vecLocations == i)
  if(length(loc) != 0){
    counter <- counter + 1
    NO3ROMSSites[counter,] <- colSums(x = NO3PC[loc,], na.rm = TRUE) / length(loc) 
    # should I average or add for nitrate?
  }
}

# NITRATE SYNCHRONY MATRIX variables:
# timesNO3 = the times that will be passed into cleandat
# matNO3 = the matrix created from the dataframe NO3ROMSSites
# cleanNO3 = the cleaned nitrate data using cleandat
# NO3Col = the number of columns in the clean nitrate data
# smNO3 = the synchrony matrix created from the clean nitrate data

# nitrate data for the synchrony matrix of nitrate data
NO3ROMSSites <- subset(NO3ROMSSites, select=-Site) # remove the 'Site' column
timesNO3 <- c(1:ncol(NO3ROMSSites)) # this is the times that will be passed into cleandat
matNO3 <- data.matrix(NO3ROMSSites) # creates the matrix from the dataframe
cleanNO3 <- cleandat(matNO3, timesNO3, clev=2)$cdat # cleans the data - 2 is enough
NO3Col <- ncol(cleanNO3)
smNO3 <- synmat(cleanNO3, 1:NO3Col, method="pearson")
diag(smNO3) <- 1
fields::image.plot(1:ncol(smNO3),1:ncol(smNO3),smNO3,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)


save.image(file = paste0(resloc, "create_synchrony.RData"))

