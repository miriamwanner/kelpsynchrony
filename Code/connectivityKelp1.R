# Miriam Wanner Connectivity and Kelp Biomass Data
# March 21, 2021

library(spatstat)
library(geosphere)

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
    kelpDataROMSSites[counter,] <- colSums(x = kelpBioPC[loc,], na.rm = TRUE) / length(loc)
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



# Probability of connectivity matrix: (1 - daily loss rate)^number of days 
# (I don't know if this is right but this is what I had written down)
dailyLossRate <- vector(,nrow(kelpDataROMSSites))
for(i in 1:nrow(kelpDataROMSSites)){
  for(j in 3:ncol(kelpDataROMSSites)){ # start at 3 because first column is the site numbers 
                                       # and take difference from adjacent points
    dailyLossRate[i] <- dailyLossRate[i] + (kelpDataROMSSites[i, j-1] - kelpDataROMSSites[i, j])
  }
}
dailyLossRate <- ((dailyLossRate/ncol(kelpDataROMSSites)) * 4) / 365 # divide by number of cols to get average
                              # and multiply by four (using quarterly data) and divide by 365 to get daily loss
probMat <- matrix(, nrow = 44, ncol = 44)
for(i in 1:length(dailyLossRate)){
  for(j in 1:nrow(oceanAvgKelpSites)){
    probMat[i, j] <- (1 - dailyLossRate[i])^oceanAvgKelpSites[i,j] # I just took this point as opposed to [j, i]
  }                                                                # I am unsure which direction to use
}






# I wanted to try to plot average dailyLossRate (of a) against dispersal time
# (if a > b then daily loss rate of a, and if a < b then daily loss rate of b)
# there are some points with daily loss rate greater than 2.4, 
# but it made it hard to see the data so I removed them (it was only about 2)
lossRateAndDispersal <- matrix(, nrow=(44*44), ncol=2)
for(i in 1:length(dailyLossRate)){
  for(j in 1:nrow(oceanAvgKelpSites)){
    if(dailyLossRate[i] < 2.4){
      lossRateAndDispersal[i, 2] = dailyLossRate[i]
      if(dailyLossRate[i] < dailyLossRate[j]){
        lossRateAndDispersal[i, 1] = oceanAvgKelpSites[i, j]
      }
    }
  }
}
plot(lossRateAndDispersal)

