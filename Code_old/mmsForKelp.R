# Miriam Wanner Connectivity and Kelp Biomass Data
# April 20, 2021

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
library(mms)
library(gtools)

# ========================================= Preparing the data ====================================================

# PREPARING KELP DATA variables:
# kelpBioOriginal = original kelp data from CSV file
# paddleCoordinates = original paddle coordinates from CSV file
# order.paddle = vector of the sites in order of paddle coordinates
# kelpBioPC = kelp data with sites ordered by paddle coordinate
# whichSitesNA = a vector listing which sites were removed (because they were all NA)

# Matrix for Kelp Biomass
kelpBioOriginal <- read.csv("/Users/miriam/Documents/Github/kelpsynchrony/Data/Kelp_Bio_2019_v3.csv")
paddleCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
order.paddle <- paddleCoordinates$Site_Number
# Biomass ordered by paddle coordinates
kelpBioPC <- kelpBioOriginal[match(order.paddle, kelpBioOriginal$Site),]
# Remove all kelp time series which consist only of NAs
# subtract 1 from the number of columns because the site column will never be NA
whichSitesNA <- which(rowSums(is.na(kelpBioPC)) == ncol(kelpBioPC) - 1)
kelpBioPC <- kelpBioPC[rowSums(is.na(kelpBioPC)) != ncol(kelpBioPC) - 1,]
# Remove all kelp time series which consist only of 0s - remove NA for the row sum
kelpBioPC <- kelpBioPC[rowSums(subset(kelpBioPC, select=-Site), na.rm = TRUE) > 0,]

# PREPARING CONNECTIVITY DATA variables:
# OceanDistanceMatrices = the matrices of connectivity data
# ocean.dist.yearly = the yearly data for the connectivity data
# ocean.dist.mainland = only the mainland data from connectivity data
# numYears = number of years the connectivity data has matrices for
# matList = the list of matrices of connectivity data for each year
# ocean.dist.avg = the average matrix of all the years of connectivity data

# Connectivity data
OceanDistanceMatrices <- readMat("/Users/miriam/Desktop/Data and code for Miriam/OceanDistanceMatrices.mat")
ocean.dist.yearly <- OceanDistanceMatrices$oceandist.yearly
ocean.dist.mainland <- ocean.dist.yearly[1:62, 1:62, ] # this is a 3D matrix
numYears <- dim(ocean.dist.mainland)[3] # this is the number of years
# for loop that adds each matrix to a list of matrices (to then average after)
matList <- list()
for(i in 1:numYears){
  matList[[i]] <- ocean.dist.mainland[,,i]
}
# from the list of matrices, take the average of each corresponding cell in each matrix for an average matrix
ocean.dist.avg <- Reduce('+', matList)/length(matList) # this is the average matrix for all time (over the time the data was collected)

# ============================================= Corresponding Kelp/ROMS sites ============================================
# Getting the kelp data and ROMS cells that are closest together/correspond to one another

# DISTANCE KELP/ROMS SITES variables:
# kelpCoordinates = coordinates of all the kelp sites from CSV file
# romsCoordinates = coordinates of all the ROMS sites from CSV file
# numKelpSites = the number of total kelp sites
# numROMSSites = the number of total ROMS sites
# lonLatKelp = two column matrix of longitude (first column) and latitude (second column) of all kelp sites
# lonLatROMS = two column matrix of longitude (first column) and latitude (second column) of all ROMS sites
# distMat = matrix of the distances between the all kelp and all ROMS sites in kilometer
# vecLocations = ROMS cell that each kelp site corresponds to (index of vector is kelp site)
# minDistForLocations = gives the distance between kelp and ROMS site (index of vector is kelp/ROMS site from vecLocations), NA if greater than 12 km away

# Get coordinates (lat and lon) of kelp and ROMS data in two data frames
kelpCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
kelpCoordinates <- kelpCoordinates[-whichSitesNA, ] # remove the sites that were removed earlier (because the entire time series was NA)
romsCoordinates <- as.data.frame(t(readMat("/Users/miriam/Desktop/Data and code for Miriam/site_centers.mat"))) # this is the center of the coordinates
# Matrix of the distance (the columns are the ROMS sites and the rows are the kelp sites)
numKelpSites <- length(as.vector(kelpCoordinates$Lat))
numROMSSites <- length(as.vector(unlist(romsCoordinates$lat))[1:62])
lonLatKelp <- matrix(c(as.vector(kelpCoordinates$Lon), as.vector(kelpCoordinates$Lat)), ncol = 2)
lonLatROMS <- matrix(c(as.vector(unlist(romsCoordinates$lon))[1:62], as.vector(unlist(romsCoordinates$lat))[1:62]), ncol = 2)
distMat <- distm(lonLatKelp, lonLatROMS, fun = distHaversine) * (0.001) # multiply to get km
vecLocations <- as.matrix(apply(distMat, 1, which.min)) # this will have the ROMS cell that each kelp site corresponds to
minDistForLocations <- vector(,length(vecLocations))
for(i in 1:length(vecLocations)){
  minDistForLocations[i] <- distMat[i, vecLocations[i]]
}
# do not use data that is greater than 12 km away (cells have 5 km radius) from the ROMS cell (make this data NA, and will remove from data set later)
for(i in 1:length(minDistForLocations)){
  if(minDistForLocations[i] > 12){
    minDistForLocations[i] <- NA
    vecLocations[i] <- NA
  }
}

# KELP/ROMS SITES CORRESPONDING TO EACHOTHER variables:
# ROMSSitesVec = vector of ROMS sites to be used, NA if no kelp data that corresponds to it
# newSites = index is the row number of the kelpDataROMSSites and value is which ROMS cell it corresponds to
# kelpDataROMSSites = kelp data (added when multiple sites fit with the same ROMS site) for the ROMS sites
# oceanAvgKelpSites = ROMS connectivity data, but deleted the sites we aren't using (no kelp site corresponds to it), diagonal is zero
# rowColToDel = the ROMS sites being deleted

# make a vector of ROMS (mainland) cells and make ones NA if there is no kelp data sites that correspond to it
# start with a vector of NA, and just change them if there is a kelp site that corresponds to it
ROMSSitesVec <- rep(NA, 62)
for(i in 1:length(vecLocations)){
  ROMSSitesVec[vecLocations[i]] <- vecLocations[i]
}
# Create the data so that the sites now correspond to one another
# new data frame that has all the kelp data, but corresponding to each of the connectivity data, and ordered by connectivity data
newSites <- vector(,43) # index is the row number of the kelpDataROMSSites and value is which ROMS cell it corresponds to 
kelpDataROMSSites <- kelpBioPC[0:43,] # first 43 rows, but the data will be replaced with ROMS data
counter <- 0
for(i in 1:nrow(kelpBioPC)){
  loc <- which(vecLocations == i)
  if(length(loc) != 0){
    counter <- counter + 1
    newSites[counter] <- i
    kelpDataROMSSites[counter,] <- colSums(x = kelpBioPC[loc,], na.rm = TRUE) # / length(loc) # only sum the elements - not average
  }
}
# delete rows/columns of the connectivity data
oceanAvgKelpSites <- ocean.dist.avg # originally is the same matrix, but then remove the sites
i <- length(ROMSSitesVec)
while(i >= 1){ # goes through the matrix and sets the row/column to all be na
  if(is.na(ROMSSitesVec[i])){
    oceanAvgKelpSites[i,] <- NA
    oceanAvgKelpSites[,i] <- NA
  }
  i <- i - 1
}
rowColToDel <- which(is.na(ROMSSitesVec))
oceanAvgKelpSites <- oceanAvgKelpSites[-rowColToDel, -rowColToDel] # removes the columns/rows that have na in them
diag(oceanAvgKelpSites) <- 0 # set the diagonal to zero

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
probMat <- (0.1)^oceanAvgKelpSites
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
# plot(x, y) # all the distances against the probability
# plot(x[x < 50], y[x < 50]) # just the distances less than 50 km against their probability (top left of previous plot)

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
# fields::image.plot(1:ncol(sm),1:ncol(sm),sm,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)

# SYMMETRIC PROBABILITY*FECUNDITY MATRIX variables:
# symProbMatF = a symmetric matrix of probMatFecundity which just takes the max of [i, j] and [j, i] from probMatFecundity
# logSym = the log10 of the symProbMatF matrix - diagonal is 1

# now I make a symmetric matrix of the probMatFecundity matrix by just using the max of [i, j] and [j, i]
symProbMatF <- pmax(probMatFecundity, t(probMatFecundity))
logSym <- log10(symProbMatF)
diag(logSym) <- 1
# fields::image.plot(1:ncol(symProbMatF),1:ncol(symProbMatF),logSym,col=heat.colors(20))
# fields::image.plot(1:ncol(logSym),1:ncol(logSym),logSym,col=heat.colors(20))

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
# fields::image.plot(1:ncol(smWaves),1:ncol(smWaves),smWaves,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)




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
# fields::image.plot(1:ncol(smNO3),1:ncol(smNO3),smNO3,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)



# USE MRM TO SEE CORRELATION BETWEEN SYNCHRONY AND SYMMETRIC PROBABILITY MATRIX:
# MRM(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(newDistMat))
MRM(as.dist(sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
MRM(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))

# MRM(as.dist(sm) ~ as.dist(logSym) + as.dist(smNO3))
# 
# MRM(as.dist(sm) ~ as.dist(newDistMat) + as.dist(logSym))
# MRM(as.dist(sm) ~ as.dist(logSym))
# MRM(as.dist(sm) ~ as.dist(newDistMat))


# USE MMS FOR MATRIX REGRESSION variables:
# matsForMR = the list of matrices used for mmsmodwts
# result = the result of using mmsmodwts on the list of matrices
# the model.names:
# 2 = distance
# 3 = connectivity (symmetric matrix for the probability matrix between sites)
# 4 = waves
# 5 = nitrate
matsForMR <- list(sm, newDistMat, symProbMatF, smWaves, smNO3)
result <- mmsmodwts(mats=matsForMR, model.names=NA, nrand=500, n=21, maxruns=50, progress=T)
result

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




# logit transform

# taking (0, 1) to (-1, 1)

diag(sm) <- NA
logit_sm <- logit(sm, min=-1, max=1)








