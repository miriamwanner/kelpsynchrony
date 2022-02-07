# Miriam Wanner Connectivity and Kelp Biomass Data
# March 31, 2021

library(spatstat)
library(geosphere)
library(R.matlab)

# Preparing the data

# Matrix for Kelp Biomass
kelpBioOriginal <- read.csv("/Users/miriam/Documents/Github/kelpsynchrony/Data/Kelp_Bio_2019_v3.csv")
paddleCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
order.paddle <- paddleCoordinates$Site_Number
# Biomass ordered by paddle coordinates
kelpBioPC <- kelpBioOriginal[match(order.paddle, kelpBioOriginal$Site),]

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





# Getting the kelp data and ROMS cells that are closest together/correspond to one another

# Get coordinates (lat and lon) of kelp and ROMS data in two data frames
kelpCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
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

# make a vector of ROMS (mainland) cells and make ones NA if there is no kelp data sites that correspond to it
# start with a vector of NA, and just change them if there is a kelp site that corresponds to it
ROMSSitesVec <- rep(NA, 62)
for(i in 1:length(vecLocations)){
  ROMSSitesVec[vecLocations[i]] <- vecLocations[i]
} # 18 sites of the ROMS data are not used
# Create the data so that the sites now correspond to one another
# new data frame that has all the kelp data, but corresponding to each of the connectivity data, and ordered by connectivity data
newSites <- vector(,44) # index is the row number of the kelpDataROMSSites and value is which ROMS cell it corresponds to 
kelpDataROMSSites <- kelpBioPC[0:44,] # first 44 rows, but the data will be replaced with ROMS data
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





# NEW CODE
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


# probability matrix - each cell is (0.1)^(days in connectivity matrix)
# set diagonal to zero
probMat <- (0.1)^oceanAvgKelpSites
diag(probMat) <- 0

# multiply each element of the probability matrix by the fecundity of each donor patch
# (similar to formula 1 in table of paper)
probMatFecundity <- probMat
for(i in 1:nrow(probMatFecundity)){
  probMatFecundity[i,] <- (probMat[i,] * kelpFecundity[i])
}


# find all distance from cells to each other in a matrix (will be symmetric)
lonLatROMS <- matrix(c(as.vector(unlist(romsCoordinates$lon))[1:62], as.vector(unlist(romsCoordinates$lat))[1:62]), ncol = 2)
newLonLatROMS <- lonLatROMS[1:44,]
count <- 1
for(i in 1:nrow(lonLatROMS)){
  if(is.element(i, newSites)){
    newLonLatROMS[count,] <- lonLatROMS[i,]
    count <- count + 1
  }
}
newDistMat <- distm(newLonLatROMS, newLonLatROMS, fun = distHaversine) * (0.001) # multiply to get km


# plot distance matrix against the probability matrix (don't plot diagonal) 
x <- as.vector(newDistMat)
y <- as.vector(log10(probMat))
plot(x[x < 50], y[x < 50])

plot(as.vector(newDistMat), as.vector(log10(probMatFecundity)))

p <- pmax(probMatFecundity, t(probMatFecundity))
a <- (probMatFecundity + t(probMatFecundity)) / 2

plot(as.vector(log10(p)), as.vector(log10(a)))

# plot p or a against the synchrony matrix

mat <- probMat
tmat <- t(probMat)

plot(as.vector(log10(mat[upper.tri(mat)])), as.vector(log10(tmat[upper.tri(tmat)])))
abline(a = 0, b = 1, col = "red")


