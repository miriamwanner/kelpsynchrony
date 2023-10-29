# script to get the new data and include the islands

setwd("/Users/miriam/Desktop/revised_kelp_code/Code") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"

# can uncomment this line if the RData has been saved:
# load(file=paste0(resloc, "new_data_processing.RData"))

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

nc_data <- nc_open('../Data/CAkelpCanopyEnv_2021.nc')

# Matrix for Kelp Biomass
kelpBioOriginal <- ncvar_get(nc_data, "biomass")
paddleCoordinates <- read.csv(file=paste0(datloc, "Paddle_coords.csv"))



# PREPARING CONNECTIVITY DATA variables:
# OceanDistanceMatrices = the matrices of connectivity data
# ocean.dist.yearly = the yearly data for the connectivity data
# ocean.dist.mainland = only the mainland data from connectivity data
# numYears = number of years the connectivity data has matrices for
# matList = the list of matrices of connectivity data for each year
# ocean.dist.avg = the average matrix of all the years of connectivity data

# Connectivity data

OceanDistanceMatrices <- readMat(con=paste0(datloc, "OceanDistanceMatrices.mat"))
ocean.dist.yearly <- OceanDistanceMatrices$oceandist.yearly
ocean.dist.mainland <- ocean.dist.yearly[1:135, 1:135, ] # this is a 3D matrix
numYears <- dim(ocean.dist.mainland)[3] # this is the number of years
# for loop that adds each matrix to a list of matrices (to then average after)
matList <- list()
for(i in 1:numYears){
  matList[[i]] <- ocean.dist.mainland[,,i]
}
# from the list of matrices, take the average of each corresponding cell in each matrix for an average matrix
ocean.dist.avg <- Reduce('+', matList)/length(matList) # this is the average matrix for all time (over the time the data was collected)

# Corresponding Kelp/ROMS sites
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
kelpLat <- ncvar_get(nc_data, "lat")
kelpLon <- ncvar_get(nc_data, "lon")
kelpSite <- 1:nrow(kelpBioOriginal)
kelpCoordinates <- data.frame(kelpLat, kelpLon, kelpSite)
names(kelpCoordinates)[1] <- "Lat"
names(kelpCoordinates)[2] <- "Lon"
names(kelpCoordinates)[3] <- "Site_Number"

romsCoordinates <- as.data.frame(t(readMat(con=paste0(datloc, "site_centers.mat")))) # this is the center of the coordinates
# Matrix of the distance (the columns are the ROMS sites and the rows are the kelp sites)
numKelpSites <- length(as.vector(kelpCoordinates$Lat))
numROMSSites <- length(as.vector(unlist(romsCoordinates$lat))[1:135])

lonLatKelp <- matrix(c(as.vector(kelpCoordinates$Lon), as.vector(kelpCoordinates$Lat)), ncol = 2)
lonLatROMS <- matrix(c(as.vector(unlist(romsCoordinates$lon))[1:135], as.vector(unlist(romsCoordinates$lat))[1:135]), ncol = 2)
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
ROMSSitesVec <- rep(NA, 135)
for(i in 1:length(vecLocations)){
  ROMSSitesVec[vecLocations[i]] <- vecLocations[i]
}
# Create the data so that the sites now correspond to one another
# new data frame that has all the kelp data, but corresponding to each of the connectivity data, and ordered by connectivity data
newSites <- vector(,120) # index is the row number of the kelpDataROMSSites and value is which ROMS cell it corresponds to 
kelpDataROMSSites <- kelpBioOriginal[0:120,] # first 120 rows, but the data will be replaced with ROMS data
counter <- 0
for(i in 1:nrow(kelpBioOriginal)){
  loc <- which(vecLocations == i)
  if(length(loc) > 1){
    counter <- counter + 1
    newSites[counter] <- i
    kelpDataROMSSites[counter,] <- colSums(x = kelpBioOriginal[loc,], na.rm = TRUE)
  }
  if(length(loc) == 1){
    counter <- counter + 1
    newSites[counter] <- i
    kelpDataROMSSites[counter,] <- kelpBioOriginal[loc,]
  }
  if (i %% 50 == 0){
    print(i)    
  }
}
# delete the row/corresponding column if there are NA values
whichSitesNA <- which(rowSums(is.na(kelpDataROMSSites)) > 0)
kelpDataROMSSites <- kelpDataROMSSites[-whichSitesNA,]

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
oceanAvgKelpSites <- oceanAvgKelpSites[-whichSitesNA, -whichSitesNA]
diag(oceanAvgKelpSites) <- 0 # set the diagonal to zero

# creating synchrony matrices

survival_rate <- 0.1 # this can be changed to recreate MRM results with different survival rates

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
kelpFecundity <- 1463 * (rowSums(sqrtKelp) / ncol(sqrtKelp))

# create a probability matrix and probability matrix multiplied by fecundity
probMat <- (survival_rate)^oceanAvgKelpSites
diag(probMat) <- 0
# multiply each element of the probability matrix by the fecundity of each donor patch
probMatFecundity <- probMat
for(i in 1:nrow(probMatFecundity)){
  probMatFecundity[i,] <- (probMat[i,] * kelpFecundity[i])
}

# DISTANCES FROM ONE CELL TO ANOTHER variables:
# lonLatROMS = matrix (two columns) of all the ROMS sites longitude (column 1) and latitude (column 2)
# newLonLatROMS = matrix (two columns) of the ROMS sites being used for kelp data longitude (column 1) and latitude (column 2)
# newDistMat = the distance matrix of the sites being used for longitude and latitude

# find all distance from cells to each other in a matrix (will be symmetric)
newLonLatROMS <- lonLatROMS[1:120,]
count <- 1
for(i in 1:nrow(lonLatROMS)){
  if(is.element(i, newSites)){
    newLonLatROMS[count,] <- lonLatROMS[i,]
    count <- count + 1
  }
}
newLonLatROMS <- newLonLatROMS[-whichSitesNA,]
newDistMat <- distm(newLonLatROMS, newLonLatROMS, fun = distHaversine) * (0.001) # multiply to get km


# SYNCHRONY MATRIX variables:
# kelpDataROMSSites = the kelp data but without the 'Site' column
# times = the times that will be passed into the cleandat function
# matKelp = the data fram kelpDataROMSSites as a matrix
# cleanKelp = the kelp data cleaned by cleandat() with a cleaning level of two
# matcol = the number of columns in cleanKelp
# sm = the synchrony matrix from using synmat()

# process the data for the synchrony matrix
times <- c(1:ncol(kelpDataROMSSites)) # this is the times that will be passed into cleandat
matKelp <- data.matrix(kelpDataROMSSites) # creates the matrix from the dataframe
cleanKelp <- cleandat(matKelp, times, clev=2)$cdat # cleans the data - 2 is enough
matcol <- ncol(cleanKelp)
sm <- synmat(cleanKelp, 1:matcol, method="pearson")
diag(sm) <- 1

# SYMMETRIC PROBABILITY*FECUNDITY MATRIX variables:
# symProbMatF = a symmetric matrix of probMatFecundity which just takes the max of [i, j] and [j, i] from probMatFecundity
# logSym = the log10 of the symProbMatF matrix - diagonal is 1

# now I make a symmetric matrix of the probMatFecundity matrix by just using the max of [i, j] and [j, i]
symProbMatF <- pmax(probMatFecundity, t(probMatFecundity))
logSym <- log10(symProbMatF)
diag(logSym) <- 1



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
waves <- ncvar_get(nc_data, "hsmax")

# get the waves sites to line up with the ROMS sites (same way as kelp sites)
wavesROMSSites <- waves[0:120,]
counter <- 0
for(i in 1:nrow(waves)){
  loc <- which(vecLocations == i)
  if(length(loc) > 1){
    counter <- counter + 1
    wavesROMSSites[counter,] <- colSums(x = waves[loc,], na.rm = TRUE) / length(loc) 
  }
  if(length(loc) == 1){
    counter <- counter + 1
    wavesROMSSites[counter,] <- waves[loc,] 
  }
  if (i %% 50 == 0){
    print(i)    
  }
}
wavesROMSSites <- wavesROMSSites[-whichSitesNA,]

# WAVES SYNCHRONY MATRIX variables:
# timesWaves = the times that will be passed into cleandat
# matWaves = the matrix created from the dataframe wavesROMSSites
# cleanWaves = the cleaned wave data using cleandat
# waveCol = the number of columns in the clean wave data
# smWaves = the synchrony matrix created from the clean wave data

# wave data for the synchrony matrix of wave data
timesWaves <- c(1:ncol(wavesROMSSites)) # this is the times that will be passed into cleandat
matWaves <- data.matrix(wavesROMSSites) # creates the matrix from the dataframe
cleanWaves <- cleandat(matWaves, timesWaves, clev=2)$cdat # cleans the data - 2 is enough
waveCol <- ncol(cleanWaves)
smWaves <- synmat(cleanWaves, 1:waveCol, method="pearson")
diag(smWaves) <- 1


# NITRATE DATA variables:
# NO3 = the original nitrate data
# NO3PC = nitrate data ordered by paddle coordinate
# NO3ROMSSites = nitrate data that matches up to the ROMS sites

# read in the nitrate data
NO3 <- ncvar_get(nc_data, "nitrate")


# get the nitrate sites to line up with the ROMS sites (same way as kelp sites)
NO3ROMSSites <- NO3[0:120,]
counter <- 0
for(i in 1:nrow(NO3)){
  loc <- which(vecLocations == i)
  if(length(loc) > 1){
    counter <- counter + 1
    NO3ROMSSites[counter,] <- colSums(x = NO3[loc,], na.rm = TRUE) / length(loc) 
  }
  if(length(loc) == 1){
    counter <- counter + 1
    NO3ROMSSites[counter,] <- NO3[loc,] 
  }
  if (i %% 50 == 0){
    print(i)    
  }
}
NO3ROMSSites <- NO3ROMSSites[-whichSitesNA,]

# NITRATE SYNCHRONY MATRIX variables:
# timesNO3 = the times that will be passed into cleandat
# matNO3 = the matrix created from the dataframe NO3ROMSSites
# cleanNO3 = the cleaned nitrate data using cleandat
# NO3Col = the number of columns in the clean nitrate data
# smNO3 = the synchrony matrix created from the clean nitrate data

# nitrate data for the synchrony matrix of nitrate data
# NO3ROMSSites <- subset(NO3ROMSSites, select=-Site) # remove the 'Site' column
timesNO3 <- c(1:ncol(NO3ROMSSites)) # this is the times that will be passed into cleandat
matNO3 <- data.matrix(NO3ROMSSites) # creates the matrix from the dataframe
cleanNO3 <- cleandat(matNO3, timesNO3, clev=2)$cdat # cleans the data - 2 is enough
NO3Col <- ncol(cleanNO3)
smNO3 <- synmat(cleanNO3, 1:NO3Col, method="pearson")
diag(smNO3) <- 1


# save for loading in other files
save.image(file=paste0(resloc, file="new_data_processing.RData"))
