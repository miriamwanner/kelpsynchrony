# script to get the new data and include the islands

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/ROMSKelpData/"
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


nc_data <- nc_open('../Data/CAkelpCanopyEnv_2021.nc')

# Matrix for Kelp Biomass
kelpBioOriginal <- ncvar_get(nc_data, "biomass")
paddleCoordinates <- read.csv(file=paste0(datloc, "Paddle_coords.csv"))
# order.paddle <- paddleCoordinates$Site_Number
# # Biomass ordered by paddle coordinates
# kelpBioPC <- kelpBioOriginal[match(order.paddle, kelpBioOriginal$Site),]
# # Remove all kelp time series which consist only of NAs
# # subtract 1 from the number of columns because the site column will never be NA
# whichSitesNA <- which(rowSums(is.na(kelpBioPC)) == ncol(kelpBioPC) - 1)
# whichSitesNA <- which(rowSums(is.na(kelpBioOriginal)) == ncol(kelpBioOriginal) - 1)
# kelpBioPC <- kelpBioPC[rowSums(is.na(kelpBioPC)) != ncol(kelpBioPC) - 1,]
# # Remove all kelp time series which consist only of 0s - remove NA for the row sum
# kelpBioPC <- kelpBioPC[rowSums(subset(kelpBioPC, select=-Site), na.rm = TRUE) > 0,]




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
# kelpCoordinates <- read.csv(file=paste0(datloc, "Paddle_coords.csv"))
kelpLat <- ncvar_get(nc_data, "lat")
kelpLon <- ncvar_get(nc_data, "lon")
kelpSite <- 1:nrow(kelpBioOriginal)
kelpCoordinates <- data.frame(kelpLat, kelpLon, kelpSite)
names(kelpCoordinates)[1] <- "Lat"
names(kelpCoordinates)[2] <- "Lon"
names(kelpCoordinates)[3] <- "Site_Number"

# kelpCoordinates <- kelpCoordinates[-whichSitesNA, ] # remove the sites that were removed earlier (because the entire time series was NA)
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
    kelpDataROMSSites[counter,] <- colSums(x = kelpBioOriginal[loc,], na.rm = TRUE) # / length(loc) # only sum the elements - not average
  }
  if(length(loc) == 1){
    counter <- counter + 1
    newSites[counter] <- i
    kelpDataROMSSites[counter,] <- kelpBioOriginal[loc,] # / length(loc) # only sum the elements - not average
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


# save as .rds or .rdata datafile
# use the command save
# when you load - the vars come into the workspace with the same names
# save rds saves one var as a file name - asign that var to a (if u want) diff var name



# for figures later
#================================================================================================
# (2) CALCULATE SEMESTERLY (SIX-MONTH) AVERAGES OF ROMS CONNECTIVITY

# Set start and end dates of analysis
t.min <- 1996.0
t.max <- 2006.5
t.range <- seq(t.min, t.max, by = 0.5)
t.range.char <- as.character(t.range)

# Recall that data represent minimum particle transit times (units are in days) 
# NOTE REMOVAL OF DIAGONAL (SELF-RECRUITMENT). This is equivalent to probability of dispersal from local cell back to itself.
OD.avg.funct <- function(ocean.dist.monthly, t.min, t.max, t.range, semester){
  sem.list <- rep(t.range, each = 6)
  sem.count <- 1:length(sem.list)
  sem.df <- data.frame(sem.list, sem.count)
  sem.range <- sem.df$sem.count[sem.df$sem.list == semester]
  OD.avg <-   (ocean.dist.monthly[,,sem.range[1]] +
                 ocean.dist.monthly[,,sem.range[2]] +
                 ocean.dist.monthly[,,sem.range[3]] +
                 ocean.dist.monthly[,,sem.range[4]] +
                 ocean.dist.monthly[,,sem.range[5]] +
                 ocean.dist.monthly[,,sem.range[6]]) / length(sem.range)
  diag(OD.avg) <- 0
  return(list(OD.avg))
}

OceanDistanceMatrices <- readMat(con=paste0(datloc, "OceanDistanceMatrices.mat"))
site_centers <- readMat(con=paste0(datloc, "site_centers.mat"))
ocean.dist.monthly <- OceanDistanceMatrices$oceandist.monthly
ocean.dist.yearly <- OceanDistanceMatrices$oceandist.yearly

OD.avgs.mats <- lapply(t.range, function(semester){
  OD.avg.funct(ocean.dist.monthly=ocean.dist.monthly, 
               t.min=t.min, 
               t.max=t.max, 
               t.range=t.range, 
               semester=semester)}
)

#================================================================================================
# (3) ORGANIZE ROMS SPATIAL DATA

# Categorize mainland and island sites
roms.cells.region <- data.frame(cell=1:135)
roms.cells.region$region <- rep(NA, 135)
roms.cells.region$region[1:62] <- c(rep("mainland", 62))
roms.cells.region$region[63:96] <- c(rep("n.islands", 34)) 
roms.cells.region$region[97:109] <- c(rep("catalina", 13))
roms.cells.region$region[110:122] <- c(rep("sanclem", 13))
roms.cells.region$region[123:130] <- c(rep("sannic", 8))
roms.cells.region$region[131:135] <- c(rep("sanbarb", 5))

# Create vectors of lat-longs for cells i and j pairs for ADJACENT PAIRS ONLY
i.roms.cell <- c(1, rep(2:61, each=2), 62, # mainland
                 rep(63:65, each=2), rep(66, times=4), rep(67:71, each=2), rep(72, times=3), rep(73:76, each=2), rep(77, times=4), rep(78, times=2), rep(79, times=3), rep(80, times=3), rep(81:86, each=2), rep(87, times=3), rep(88, times=2), rep(89, time=3), rep(90, times=3), rep(91, times=2), rep(92, times=3), rep(93, times=3), rep(94:96, each=2),      #north islands
                 rep(97:109, each=2), # catalina
                 rep(110:122, each=2), # sanclem
                 rep(123:130, each=2), # sannic
                 rep(131:135, each=2) # sanbarb
) 

j.roms.cell <- c(c(rbind(2:62,1:61)),96,64,63,65,64,66,65,67,92,93,66,68,67,69,68,70,69,71,70,72,71,73,88,72,74,73,75,74,76,75,77,76,78,79,80,77,79,77,78,80,77,79,81,80,82,81,83,82,84,83,85,84,86,85,87,86,89,90,72,89,87,88,90,87,89,91,90,92,66,91,93,66,92,94,93,95,94,96,63,95,c(109,rbind(98:109,97:108),97),c(122,rbind(111:122,110:121),110),c(130,rbind(124:130,123:129),123),c(135,rbind(132:135,131:134),131))  

roms.pairs <- data.frame(i.roms.cell, j.roms.cell)
roms.pairs$i.roms.lat <- mapvalues(roms.pairs$i.roms.cell, from=1:135, to=c(t(site_centers$lat)))
roms.pairs$i.roms.long <- mapvalues(roms.pairs$i.roms.cell, from=1:135, to=c(t(site_centers$lon)))
roms.pairs$j.roms.lat <- mapvalues(roms.pairs$j.roms.cell, from=1:135, to=c(t(site_centers$lat)))
roms.pairs$j.roms.long <- mapvalues(roms.pairs$j.roms.cell, from=1:135, to=c(t(site_centers$lon)))
roms.pairs$region <- mapvalues(roms.pairs$i.roms.cell, from=1:135, to=roms.cells.region$region)

# Measure "great circle" distance between adjacent mainland sites (in km)
i.long.lat.mat <- as.matrix(cbind(roms.pairs$i.roms.long, roms.pairs$i.roms.lat))
j.long.lat.mat <- as.matrix(cbind(roms.pairs$j.roms.long, roms.pairs$j.roms.lat))

dist.funct <- function(i.long.lat.mat, j.long.lat.mat, r){distm(x=i.long.lat.mat[r,], y=j.long.lat.mat[r,])/1000}

roms.pairs$dist <- sapply(1:dim(roms.pairs)[1], function(r){dist.funct(i.long.lat.mat=i.long.lat.mat, j.long.lat.mat=j.long.lat.mat, r=r)})

roms.pairs <- roms.pairs[,c(7,1:6,8)] # Rearrange

# First, make a function: For a given i -> j pair of cells, return the number of days
ij.pairs.list <- cbind(roms.pairs$i.roms.cell, roms.pairs$j.roms.cell)
ij.days <- function(OD.mat, ij.pairs){as.matrix(as.data.frame(OD.mat))[[ij.pairs[1], ij.pairs[2]]]}

# For each semester, count the number of days between all i -> j cells 
all.ij.days <- function(ij.pairs.list, OD.mat){
  apply(X=ij.pairs.list, MARGIN=1, FUN=function(row){ij.days(OD.mat=OD.mat, ij.pairs=row)})
}

# Duplicate data frame for all semesters (22)
roms.pairs.list <- list(roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs, roms.pairs)

# Use lapply to run all.ij.days function across all semesters
days.list <- lapply(1:22, function(l){all.ij.days(ij.pairs.list=ij.pairs.list, OD.mat=OD.avgs.mats[l])})

# Merge 2 lists together into roms.pairs.list, and divide to get days per km
for(x in 1:22){
  roms.pairs.list[[x]]$days <- days.list[[x]]
  roms.pairs.list[[x]]$days.per.km <- roms.pairs.list[[x]]$days / roms.pairs.list[[x]]$dist
}

save(roms.pairs.list, file=here::here("roms.pairs.list.rda"))

# Split regions into separate lists of data frames
roms.pairs.list.mainland <- lapply(roms.pairs.list, subset, region=="mainland")
roms.pairs.list.n.islands <- lapply(roms.pairs.list, subset, region=="n.islands")
roms.pairs.list.catalina <- lapply(roms.pairs.list, subset, region=="catalina")
roms.pairs.list.sanclem <- lapply(roms.pairs.list, subset, region=="sanclem")
roms.pairs.list.sannic <- lapply(roms.pairs.list, subset, region=="sannic")
roms.pairs.list.sanbarb <- lapply(roms.pairs.list, subset, region=="sanbarb")

# Split again based on north vs south direction, or clockwise vs counterclockwise direction
roms.pairs.list.mainland.SN <- lapply(roms.pairs.list.mainland, subset, i.roms.cell < j.roms.cell)
roms.pairs.list.mainland.NS <- lapply(roms.pairs.list.mainland, subset, i.roms.cell > j.roms.cell)

# Measure cumulative distance from southernmost cell (in km)
for(x in 1:22){
  roms.pairs.list.mainland.SN[[x]]$cum.dist <- cumsum(roms.pairs.list.mainland.SN[[x]]$dist)
  roms.pairs.list.mainland.NS[[x]]$cum.dist <- cumsum(roms.pairs.list.mainland.SN[[x]]$dist)
}

# Merge list of data frames together
ALL.roms.pairs.list.mainland.SN <- do.call("rbind", roms.pairs.list.mainland.SN)
ALL.roms.pairs.list.mainland.SN$year <- rep(t.range, each=length(unique(ALL.roms.pairs.list.mainland.SN$i.roms.cell)))
ALL.roms.pairs.list.mainland.NS <- do.call("rbind", roms.pairs.list.mainland.NS)
ALL.roms.pairs.list.mainland.NS$year <- ALL.roms.pairs.list.mainland.SN$year 

roms.pairs.list.mainland.SN[[1]]$long <- rowMeans(cbind(roms.pairs.list.mainland.SN[[1]]$i.roms.long, roms.pairs.list.mainland.SN[[1]]$j.roms.long))
roms.pairs.list.mainland.SN[[1]]$lat <- rowMeans(cbind(roms.pairs.list.mainland.SN[[1]]$i.roms.lat, roms.pairs.list.mainland.SN[[1]]$j.roms.lat))
roms.pairs.list.mainland.NS[[1]]$long <- rowMeans(cbind(roms.pairs.list.mainland.NS[[1]]$i.roms.long, roms.pairs.list.mainland.SN[[1]]$j.roms.long))
roms.pairs.list.mainland.NS[[1]]$lat <- rowMeans(cbind(roms.pairs.list.mainland.NS[[1]]$i.roms.lat, roms.pairs.list.mainland.SN[[1]]$j.roms.lat))

# roms.pairs.list.subset <- subset(roms.pairs.list.mainland.SN[[1]], !is.na(ROMSSitesVec))
# roms.pairs.list.subset <- subset()
# roms.pairs.list.subset <- subset(roms.pairs.list[[1]], !is.na(ROMSSitesVec))
roms.pairs.list.subset <- data.frame(kelpClust$coords$lon, kelpClust$coords$lat, kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"

# save.image(file = paste0(resloc, "get_ROMS.RData"))

# ========================================= creating synchrony matrices =========================================================================

survival_rate <- 0.1

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
# lonLatROMS <- matrix(c(as.vector(unlist(romsCoordinates$lon))[1:62], as.vector(unlist(romsCoordinates$lat))[1:62]), ncol = 2)
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
# kelpDataROMSSites <- subset(kelpDataROMSSites, select=-Site) # remove the 'Site' column
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
# waves <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Hs_Max_2019.csv")
# wavesPC <- waves[match(order.paddle, waves$Site),]
waves <- ncvar_get(nc_data, "hsmax")

# get the waves sites to line up with the ROMS sites (same way as kelp sites)
wavesROMSSites <- waves[0:120,]
counter <- 0
for(i in 1:nrow(waves)){
  loc <- which(vecLocations == i)
  if(length(loc) > 1){
    counter <- counter + 1
    wavesROMSSites[counter,] <- colSums(x = waves[loc,], na.rm = TRUE) / length(loc) 
    # should I average or add for waves? I think average because it is with height
  }
  if(length(loc) == 1){
    counter <- counter + 1
    wavesROMSSites[counter,] <- waves[loc,] 
    # should I average or add for waves? I think average because it is with height
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
# wavesROMSSites <- subset(wavesROMSSites, select=-Site) # remove the 'Site' column
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
# NO3 <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/NO3_Mean_2019.csv")
# NO3PC <- NO3[match(order.paddle, NO3$Site),]
NO3 <- ncvar_get(nc_data, "nitrate")


# get the nitrate sites to line up with the ROMS sites (same way as kelp sites)
NO3ROMSSites <- NO3[0:120,]
counter <- 0
for(i in 1:nrow(NO3)){
  loc <- which(vecLocations == i)
  if(length(loc) > 1){
    counter <- counter + 1
    NO3ROMSSites[counter,] <- colSums(x = NO3[loc,], na.rm = TRUE) / length(loc) 
    # should I average or add for nitrate?
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
fields::image.plot(1:ncol(smNO3),1:ncol(smNO3),smNO3,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)


# save.image(file = paste0(resloc, "create_synchrony.RData"))

# ===================================== figures ===============================================

# 8 models with MRM

diag(sm) <- NA
logit_sm <- logit(sm, min=-1, max=1)
diag(sm) <- 1

source("altered_mrm_function.R")

# transport, linear, logit
t_lin_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, linear, not logit
t_lin_nlogit <- my_mrm(as.dist(sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, log, logit
t_log_logit <- MRM(as.dist(logit_sm) ~ as.dist(logSymProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, log, not logit
t_log_nlogit <- MRM(as.dist(sm) ~ as.dist(logSymProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, linear, logit
c_lin_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, linear, not logit
c_lin_nlogit <- my_mrm(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, log, logit
c_log_logit <- MRM(as.dist(logit_sm) ~ as.dist(logSym) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, log, not logit
c_log_nlogit <- MRM(as.dist(sm) ~ as.dist(logSym) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))

predictors <- c("dispersal", "distance", "waves", "nitrate")

all_sr_0.02 <- data.frame(predictors,
                     t_lin_logit$coef[7:10],
                     t_lin_nlogit$coef[7:10],
                     t_log_logit$coef[7:10],
                     t_log_nlogit$coef[7:10],
                     c_lin_logit$coef[7:10],
                     c_lin_nlogit$coef[7:10],
                     c_log_logit$coef[7:10],
                     c_log_nlogit$coef[7:10])
saveRDS(all_sr_0.02, file="all_sr_0.02.Rds")


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
# arrows(75, 0.6, 25, 0.4, length=0.1, angle=20, lwd=2)
# text(100, 0.7, "fast decline for distances")
# text(100, 0.64, "less than ~50km")
# arrows(300, 0.3, 250, 0.15, length=0.1, angle=20, lwd=2)
# text(300, 0.4, "slow decline for distances")
# text(300, 0.34, "greater than ~50km")
dev.off()

# plots corrected for other variables


pdf(file=paste0(resloc, "transparent_plots_figure.pdf"), width=7, height=4.2)
par(mfrow=c(2, 3), mai = c(0.4, 0.4, 0.4, 0.1), oma=c(1.5, 1.5, 0, 0))
# layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))

transport_vec <- symProbMat[upper.tri(symProbMat, diag=FALSE)]
log_transport_vec <- logSymProbMat[upper.tri(logSymProbMat, diag=FALSE)]
# plot(transport_vec, synch_vec, main="Transport", xlab="Transport", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.2))
plot(log_transport_vec, synch_vec, main="Dispersal (log scale)", xlab="log(Transport)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))

connectivity_vec <- symProbMatF[upper.tri(symProbMatF, diag=FALSE)]
log_connectivity_vec <- logSym[upper.tri(logSym, diag=FALSE)]
# plot(connectivity_vec, synch_vec, main="Connectivity", xlab="Connectivity", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.2))
plot(log_connectivity_vec, synch_vec, main="Connectivity (log scale)", xlab="log(Connectivity)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))

dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)]
plot(dist_vec, synch_vec, main="Distance (km)", xlab="Distance (km)", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))

waves_vec <- smWaves[upper.tri(smWaves, diag=FALSE)]
plot(waves_vec, synch_vec, main="Synchrony in wave height", xlab="Waves", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))

no3_vec <- smNO3[upper.tri(smNO3, diag=FALSE)]
plot(no3_vec, synch_vec, main="Synchrony in nitrate", xlab="Nitrate", ylab="Synchrony", pch=20, col=rgb(red = 0, green = 0, blue = 0, alpha=0.04))

mtext(outer=TRUE, "Predictor", side=1)
mtext(outer=TRUE, "Synchrony", side=2)
dev.off()



# plots (same as above), but using ghex
conn_df <- data.frame(connectivity_vec, synch_vec)
log_conn_df <- data.frame(log_connectivity_vec, synch_vec)
ggplot(conn_df, aes(connectivity_vec, synch_vec)) + xlab("Connectivity") + ylab("Synchrony") + geom_hex()
ggplot(log_conn_df, aes(log_connectivity_vec, synch_vec)) + xlab("log(Connectivity)") + ylab("Synchrony") + geom_hex()

transport_df <- data.frame(transport_vec, synch_vec)
log_transport_df <- data.frame(log_transport_vec, synch_vec)
ggplot(transport_df, aes(transport_vec, synch_vec)) + xlab("Transport") + ylab("Synchrony") + geom_hex()
ggplot(log_transport_df, aes(log_transport_vec, synch_vec)) + xlab("log(Transport)") + ylab("Synchrony") + geom_hex()

waves_df <- data.frame(waves_vec, synch_vec)
ggplot(waves_df, aes(waves_vec, synch_vec)) + xlab("Waves") + ylab("Synchrony") + geom_hex()

no3_df <- data.frame(no3_vec, synch_vec)
ggplot(no3_df, aes(no3_vec, synch_vec)) + xlab("Nitrate") + ylab("Synchrony") + geom_hex()

dist_df <- data.frame(dist_vec, synch_vec)
ggplot(dist_df, aes(dist_vec, synch_vec)) + xlab("Distance") + ylab("Synchrony") + geom_hex()



# clustering map

times <- c(1:ncol(kelpDataROMSSites))

colnames(newLonLatROMS)[2] <- "lat"
colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
coords <- as.data.frame(newLonLatROMS)

cleanKelp <- cleandat(data.matrix(kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp)
kelpClust <- clust(dat=cleanKelp,times=1:matcol,coords=coords,method="pearson")
cluster_numbers <- get_clusters(kelpClust)[[2]]
pdf(file=paste0(resloc, "new_clustering.pdf"), width=5, height=4)
plotmap(kelpClust)
dev.off()

# code sent from max
# ---------------------------------------------------------------------------------------------------
# Replace these with wherever you store the coastline data
coastline.data.dir        <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_f.b")
coastline.border.data.dir <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_f.b")
inset.data.dir            <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_l.b")
inset.border.data.dir     <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_l.b")

# ---------------------------------------------------------------------------------------------------
# Import coastline data from NOAA
coast.limits <- list(x = c(-125 + 360, -114 + 360), y = c(30, 39))

# Coastline polygons
coast.polys <- importGSHHS(gshhsDB=coastline.data.dir,  
                           xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.polys$X <- coast.polys$X - 360  # Get longitude back into units of degrees east
coast.polys.sp <- maptools::PolySet2SpatialPolygons(coast.polys, close_polys=FALSE) # Convert to spatial polygons
polys <- fortify(coast.polys)

# Coastline borders
coast.border <- importGSHHS(gshhsDB=coastline.border.data.dir,
                            xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.border$X <- coast.border$X - 360  # Get longitude back into units of degrees east
coast.border.sp <- maptools::PolySet2SpatialPolygons(coast.border, close_polys=FALSE) # Convert to spatial polygons
border <- fortify(coast.border)

# ---------------------------------------------------------------------------------------------------
# Set parameters for maps
lat.0 <- 32.45
lat.1 <- 34.65

long.0 <- -120.65
long.1 <- -117

# ---------------------------------------------------------------------------------------------------
# Create scaled-out inset for map
inset.limits = list(x = c(-180 + 360, -60 + 360), y = c(0, 90))
inset.polys = importGSHHS(gshhsDB=inset.data.dir, xlim=inset.limits$x, ylim=inset.limits$y, maxLevel=1, n=0)
inset.polys = fortify(inset.polys) 

inset.border = importGSHHS(gshhsDB=inset.border.data.dir, xlim=inset.limits$x, ylim=inset.limits$y, maxLevel=1, n=0)
inset.border = fortify(inset.border)

inset.map <- ggplot() + 
  coord_map(xlim=c(-121,-72.5), ylim=c(23,50), projection="lambert", lat0=23, lat1=50) +
  geom_polygon(data=inset.polys, aes(x=X, y=Y, group=PID), fill="burlywood3") +
  geom_path(data=inset.polys, aes(x=X, y=Y, group=PID), color="burlywood4", size=0.1) +
  #geom_path(data=inset.border, aes(x=X, y=Y, group=PID), color="black", size=0.5) +
  annotate(geom = 'rect', ymin=lat.0, ymax=lat.1, xmin=long.0, xmax=long.1, color="black", fill="transparent", size=0.7) +
  theme_bw() + 
  labs(x=NULL,y=NULL) + 
  theme(
    #panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background =  element_rect(fill="lightblue1"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA, size = 1),
    plot.margin = rep(unit(0,"null"),4),
    panel.margin = unit(0,"null"),
    axis.ticks.length = unit(0,"null"),
    axis.ticks.margin = unit(0,"null"),
    text = element_text(size=22),
    axis.ticks = element_blank(), 
    axis.text.x = element_blank(), axis.title.x =element_blank(), 
    axis.title.y= element_blank(), axis.text.y = element_blank()) 
inset.map

# ---------------------------------------------------------------------------------------------------
base.map <-
  ggplot() +
  coord_map(xlim=c(long.0, long.1), ylim=c(lat.0, lat.1), projection="lambert", lat0=lat.0, lat1=lat.1) +
  geom_polygon(data=polys, aes(x=X, y=Y, group=PID), fill="navajowhite2") +
  geom_path(data=polys, aes(x=X, y=Y, group=PID), color="navajowhite4", size=0.1) +
  geom_path(data=border, aes(x=X, y=Y, group=PID), color="navajowhite4", size=0.1) +
  ggsn::scalebar(x.min = long.0, x.max = long.1, y.min = lat.0, y.max = lat.1,  
                 transform = TRUE, model = 'WGS84',
                 location = "bottomleft", dist = 50, dist_unit = "km",
                 st.bottom = FALSE, height=0.02, st.dist = 0.035, st.size = 4,
                 st.color = 'black', box.fill = 'black', box.color = 'transparent') +
  scale_x_continuous(breaks=rev(seq(floor(long.0), ceiling(long.1), by = 1)), 
                     labels= -1 * rev(seq(floor(long.0), ceiling(long.1), by = 1))) +
  scale_y_continuous(breaks= seq(floor(lat.0), ceiling(lat.1), by=1),
                     labels= seq(floor(lat.0), ceiling(lat.1), by = 1)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +
  ylab(expression(paste(Latitude," (",degree,N,")", sep=""))) +
  annotate(geom = 'text', x = -119.75, y = 32.9, label = 'Pacific Ocean', size = 6) +
  annotate(geom = 'text', x = -117.7, y = 34.6, label = 'California, USA', size = 6)+
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightblue1"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1.5),
        panel.spacing = unit(1,"lines"),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        text = element_text(size = 20),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 20),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 20),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size=1, color="black"),
        axis.ticks.length = unit(0.25, "cm")) 
base.map

# ---------------------------------------------------------------------------------------------------
# Add polygons of ROMS cells to base map

# MODIFY THIS BASED ON YOUR EXISTING CODE
roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=9, shape=21, color="grey") +
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1, color="black")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(axis.text.x = element_text(color="black")) +
  theme(axis.text.y = element_text(color="black")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2))
roms.map


# getting the different colors of clusters
# for mainland data:
# roms.pairs.list.subset$Region <- rep(NA, 117)
# roms.pairs.list.subset$Region[which(cluster_numbers == 1)] <- c(rep("Northerly", length(which(cluster_numbers == 1))))
# roms.pairs.list.subset$Region[which(cluster_numbers == 2)] <- c(rep("Central", length(which(cluster_numbers == 2))))
# roms.pairs.list.subset$Region[which(cluster_numbers == 3)] <- c(rep("Southerly", length(which(cluster_numbers == 3))))
# for all data:
# roms.pairs.list.subset$region
# roms.pairs.list.subset$region[which(cluster_numbers == 1)] <- c(rep("Cluster 1", length(which(cluster_numbers == 1))))
# roms.pairs.list.subset$region[which(cluster_numbers == 2)] <- c(rep("Cluster 2", length(which(cluster_numbers == 2))))
# roms.pairs.list.subset$i.roms.lat <- newLonLatROMS[,2]
# roms.pairs.list.subset$i.roms.long <- newLonLatROMS[,1]

roms.pairs.list.subset <- data.frame(kelpClust$coords$lon, kelpClust$coords$lat, kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "all_clustering_map.pdf"), width=7, height=5)
roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
  # geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  # geom_point(data=kelpCoordinates, aes(x=Lon, y=Lat), size=0.5, color="grey") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1, color="black")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(axis.text.x = element_text(color="black")) +
  theme(axis.text.y = element_text(color="black")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2))
roms.map
dev.off()



# colorful matrices


pdf(file=paste0(resloc, "new_matrices_figure.pdf"), width=7, height=10)
par(mfrow=c(2, 3), mai = c(0.6, 0.6, 0.6, 0.8), oma=c(2, 4, 0, 0))
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
n <- 100
# colorful matrix of sm
fields::image.plot(1:ncol(sm),1:ncol(sm),sm,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(a) Kelp Synchrony",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of smWaves
fields::image.plot(1:ncol(smWaves),1:ncol(smWaves),smWaves,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(b) Wave Synchrony",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of smNO3
fields::image.plot(1:ncol(smNO3),1:ncol(smNO3),smNO3,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(c) Nitrate Synchrony",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of newDistMat
fields::image.plot(1:ncol(newDistMat),1:ncol(newDistMat),newDistMat,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(d) Distance (km)",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of transport - log logSymProbMat
diag(logSymProbMat) <- 0
fields::image.plot(1:ncol(logSymProbMat),1:ncol(logSymProbMat),logSymProbMat,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(e) Log(Transport Index)",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of connectivity - log logSym
fields::image.plot(1:ncol(logSym),1:ncol(logSym),logSym,col=hcl.colors(n, palette="viridis"),xlab="",ylab="",main="(f) Log(Connectivity Index)",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
mtext(outer=TRUE, "Site Index", side=1)
mtext(outer=TRUE, "Site Index", side=2)
dev.off()


save.image(file = paste0(resloc, "sr0.5_new_data_processing.RData"))





# MRM for only islands (sites 50-117)

sm_islands <- sm[-(1:50), -(1:50)]
symProbMat_islands <- symProbMat[-(1:50), -(1:50)]
logSymProbMat_islands <- logSymProbMat[-(1:50), -(1:50)]
symProbMatF_islands <- symProbMatF[-(1:50), -(1:50)]
logSym_islands <- logSym[-(1:50), -(1:50)]
newDistMat_islands <- newDistMat[-(1:50), -(1:50)]
smWaves_islands <- smWaves[-(1:50), -(1:50)]
smNO3_islands <- smNO3[-(1:50), -(1:50)]


diag(sm_islands) <- NA
logit_sm_islands <- logit(sm_islands, min=-1, max=1)
diag(sm_islands) <- 1

source("altered_mrm_function.R")

# transport, linear, logit
t_lin_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(symProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, linear, not logit
t_lin_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(symProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, log, logit
t_log_logit <- MRM(as.dist(logit_sm_islands) ~ as.dist(logSymProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, log, not logit
t_log_nlogit <- MRM(as.dist(sm_islands) ~ as.dist(logSymProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, linear, logit
c_lin_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(symProbMatF_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, linear, not logit
c_lin_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(symProbMatF_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, log, logit
c_log_logit <- MRM(as.dist(logit_sm_islands) ~ as.dist(logSym_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, log, not logit
c_log_nlogit <- MRM(as.dist(sm_islands) ~ as.dist(logSym_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))

predictors <- c("dispersal", "distance", "waves", "nitrate")

islands_sr_0.5 <- data.frame(predictors,
                          t_lin_logit$coef[7:10],
                          t_lin_nlogit$coef[7:10],
                          t_log_logit$coef[7:10],
                          t_log_nlogit$coef[7:10],
                          c_lin_logit$coef[7:10],
                          c_lin_nlogit$coef[7:10],
                          c_log_logit$coef[7:10],
                          c_log_nlogit$coef[7:10])
saveRDS(islands_sr_0.5, file="islands_sr_0.5.Rds")






# islands clustering map

island_kelpDataROMSSites <- kelpDataROMSSites[51:117,]

times <- c(1:ncol(island_kelpDataROMSSites))

colnames(newLonLatROMS)[2] <- "lat"
colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
island_coords <- as.data.frame(newLonLatROMS[51:117,])

island_cleanKelp <- cleandat(data.matrix(island_kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
island_matcol <- ncol(island_cleanKelp)
island_kelpClust <- clust(dat=island_cleanKelp,times=1:island_matcol,coords=island_coords,method="pearson")
island_cluster_numbers <- get_clusters(island_kelpClust)[[2]]
pdf(file=paste0(resloc, "island_clustering.pdf"), width=5, height=4)
plotmap(island_kelpClust)
dev.off()

# code sent from max
# ---------------------------------------------------------------------------------------------------
# Replace these with wherever you store the coastline data
coastline.data.dir        <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_f.b")
coastline.border.data.dir <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_f.b")
inset.data.dir            <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_l.b")
inset.border.data.dir     <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_l.b")

# ---------------------------------------------------------------------------------------------------
# Import coastline data from NOAA
coast.limits <- list(x = c(-125 + 360, -114 + 360), y = c(30, 39))

# Coastline polygons
coast.polys <- importGSHHS(gshhsDB=coastline.data.dir,  
                           xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.polys$X <- coast.polys$X - 360  # Get longitude back into units of degrees east
coast.polys.sp <- maptools::PolySet2SpatialPolygons(coast.polys, close_polys=FALSE) # Convert to spatial polygons
polys <- fortify(coast.polys)

# Coastline borders
coast.border <- importGSHHS(gshhsDB=coastline.border.data.dir,
                            xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.border$X <- coast.border$X - 360  # Get longitude back into units of degrees east
coast.border.sp <- maptools::PolySet2SpatialPolygons(coast.border, close_polys=FALSE) # Convert to spatial polygons
border <- fortify(coast.border)

# ---------------------------------------------------------------------------------------------------
# Set parameters for maps
lat.0 <- 32.45
lat.1 <- 34.65

long.0 <- -120.65
long.1 <- -117

# ---------------------------------------------------------------------------------------------------
# Create scaled-out inset for map
inset.limits = list(x = c(-180 + 360, -60 + 360), y = c(0, 90))
inset.polys = importGSHHS(gshhsDB=inset.data.dir, xlim=inset.limits$x, ylim=inset.limits$y, maxLevel=1, n=0)
inset.polys = fortify(inset.polys) 

inset.border = importGSHHS(gshhsDB=inset.border.data.dir, xlim=inset.limits$x, ylim=inset.limits$y, maxLevel=1, n=0)
inset.border = fortify(inset.border)

inset.map <- ggplot() + 
  coord_map(xlim=c(-121,-72.5), ylim=c(23,50), projection="lambert", lat0=23, lat1=50) +
  geom_polygon(data=inset.polys, aes(x=X, y=Y, group=PID), fill="burlywood3") +
  geom_path(data=inset.polys, aes(x=X, y=Y, group=PID), color="burlywood4", size=0.1) +
  #geom_path(data=inset.border, aes(x=X, y=Y, group=PID), color="black", size=0.5) +
  annotate(geom = 'rect', ymin=lat.0, ymax=lat.1, xmin=long.0, xmax=long.1, color="black", fill="transparent", size=0.7) +
  theme_bw() + 
  labs(x=NULL,y=NULL) + 
  theme(
    #panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background =  element_rect(fill="lightblue1"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA, size = 1),
    plot.margin = rep(unit(0,"null"),4),
    panel.margin = unit(0,"null"),
    axis.ticks.length = unit(0,"null"),
    axis.ticks.margin = unit(0,"null"),
    text = element_text(size=22),
    axis.ticks = element_blank(), 
    axis.text.x = element_blank(), axis.title.x =element_blank(), 
    axis.title.y= element_blank(), axis.text.y = element_blank()) 
inset.map

# ---------------------------------------------------------------------------------------------------
base.map <-
  ggplot() +
  coord_map(xlim=c(long.0, long.1), ylim=c(lat.0, lat.1), projection="lambert", lat0=lat.0, lat1=lat.1) +
  geom_polygon(data=polys, aes(x=X, y=Y, group=PID), fill="navajowhite2") +
  geom_path(data=polys, aes(x=X, y=Y, group=PID), color="navajowhite4", size=0.1) +
  geom_path(data=border, aes(x=X, y=Y, group=PID), color="navajowhite4", size=0.1) +
  ggsn::scalebar(x.min = long.0, x.max = long.1, y.min = lat.0, y.max = lat.1,  
                 transform = TRUE, model = 'WGS84',
                 location = "bottomleft", dist = 50, dist_unit = "km",
                 st.bottom = FALSE, height=0.02, st.dist = 0.035, st.size = 4,
                 st.color = 'black', box.fill = 'black', box.color = 'transparent') +
  scale_x_continuous(breaks=rev(seq(floor(long.0), ceiling(long.1), by = 1)), 
                     labels= -1 * rev(seq(floor(long.0), ceiling(long.1), by = 1))) +
  scale_y_continuous(breaks= seq(floor(lat.0), ceiling(lat.1), by=1),
                     labels= seq(floor(lat.0), ceiling(lat.1), by = 1)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +
  ylab(expression(paste(Latitude," (",degree,N,")", sep=""))) +
  annotate(geom = 'text', x = -119.75, y = 32.9, label = 'Pacific Ocean', size = 6) +
  annotate(geom = 'text', x = -117.7, y = 34.6, label = 'California, USA', size = 6)+
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightblue1"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1.5),
        panel.spacing = unit(1,"lines"),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        text = element_text(size = 20),
        axis.text.x = element_text(margin=margin(0.25,0.25,0.5,0.25, "cm"), color = "black", size = 20),
        axis.text.y = element_text(margin=margin(0.25,0.25,0.25,0.25, "cm"), color = "black", size = 20),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        axis.ticks = element_line(size=1, color="black"),
        axis.ticks.length = unit(0.25, "cm")) 
base.map

# ---------------------------------------------------------------------------------------------------
# Add polygons of ROMS cells to base map

# MODIFY THIS BASED ON YOUR EXISTING CODE
roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=9, shape=21, color="grey") +
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1, color="black")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(axis.text.x = element_text(color="black")) +
  theme(axis.text.y = element_text(color="black")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2))
roms.map


# roms.pairs.list.subset$Cluster <- c(rep(NA, length(roms.pairs.list.subset$region)))
# roms.pairs.list.subset$Cluster[which(island_cluster_numbers == 1) + 50] <- c(rep("Cluster 1", length(which(island_cluster_numbers == 1))))
# roms.pairs.list.subset$Cluster[which(island_cluster_numbers == 2) + 50] <- c(rep("Cluster 2", length(which(island_cluster_numbers == 2))))

roms.pairs.list.subset <- data.frame(island_kelpClust$coords$lon, island_kelpClust$coords$lat, island_kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "island_clustering_map.pdf"), width=7, height=5)
roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
  # geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  # geom_point(data=kelpCoordinates, aes(x=Lon, y=Lat), size=0.5, color="grey") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1, color="black")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(axis.text.x = element_text(color="black")) +
  theme(axis.text.y = element_text(color="black")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2))
roms.map
dev.off()




# mainland clustering map
main_kelpDataROMSSites <- kelpDataROMSSites[1:49,]

times <- c(1:ncol(main_kelpDataROMSSites))

colnames(newLonLatROMS)[2] <- "lat"
colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
main_coords <- as.data.frame(newLonLatROMS[1:49,])

main_cleanKelp <- cleandat(data.matrix(main_kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
main_matcol <- ncol(main_cleanKelp)
main_kelpClust <- clust(dat=main_cleanKelp,times=1:main_matcol,coords=main_coords,method="pearson")
main_cluster_numbers <- get_clusters(main_kelpClust)[[2]]
# recode(main_cluster_numbers, `1`=2, `2`=1)
pdf(file=paste0(resloc, "main_clustering.pdf"), width=5, height=4)
plotmap(main_kelpClust)
dev.off()

roms.pairs.list.subset <- data.frame(main_kelpClust$coords$lon, main_kelpClust$coords$lat, main_kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "new_mainland_clustering_map.pdf"), width=7, height=5)
roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
  # geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  # geom_point(data=kelpCoordinates, aes(x=Lon, y=Lat), size=0.5, color="grey") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
  theme(legend.title = element_text(size=18)) +
  theme(axis.ticks = element_line(size=1, color="black")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.ticks.margin = unit(0.25, "cm")) +
  theme(axis.text.x = element_text(color="black")) +
  theme(axis.text.y = element_text(color="black")) +
  theme(plot.title = element_text(size=25)) +
  theme(text = element_text(size=22)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +  
  ylab(expression(paste(Latitude," (",degree,N,")", sep="")))  +  
  theme(axis.title.x=element_text(vjust=-0.7)) +  
  theme(axis.title.y=element_text(vjust=1.2))
roms.map
dev.off()




# average time series plots
# colors cluster 1: "#F8766D", cluster 2: "#00BFC4"

# years 
# years_kelp <- ncvar_get(nc_data, "year")


# for all data
all_clust1 <- which(cluster_numbers == 1)
all_clust2 <- which(cluster_numbers == 2)

all_clust1_TS <- cleanKelp[-c(all_clust2),]
all_clust2_TS <- cleanKelp[-c(all_clust1),]

all_clust1_avg <- colMeans(all_clust1_TS)
all_clust2_avg <- colMeans(all_clust2_TS)

# all_x <- 1:length(all_clust1_avg)
all_x <- (1:length(all_clust1_avg) * 0.25) + 1983.75
# plot(main_x, mainland_clust1_avg, xlim=range(main_x), 
#      ylim=range(mainland_clust1_avg), xlab="Time", ylab="Kelp Biomass",
#      main="Mainland Average Time Series", pch=16, cex.lab=1.5, cex.main=1.8, cex.axis=1.8)

pdf(file=paste0(resloc, "all_avg_ts.pdf"), width=7, height=5)
plot(1, type="n", xlim=range(all_x), 
     ylim=range(all_clust2_avg), xlab="Year", ylab="Kelp biomass",
     # main="Time series of average kelp biomass: all locations", 
     cex.lab=1.7, cex.main=1.8, cex.axis=1.7,
     font.main=1)
lines(all_x[order(all_x)], all_clust1_avg[order(all_x)],
      xlim=range(all_x), ylim=range(all_clust1_avg), pch=16, col="#F8766D", lwd="2")
lines(all_x[order(all_x)], all_clust2_avg[order(all_x)],
      xlim=range(all_x), ylim=range(all_clust2_avg), pch=16, col="#00BFC4", lwd="2")
dev.off()

# for island data
island_clust1 <- which(island_cluster_numbers == 1)
island_clust2 <- which(island_cluster_numbers == 2)

island_clust1_TS <- island_cleanKelp[-c(island_clust2),]
island_clust2_TS <- island_cleanKelp[-c(island_clust1),]

island_clust1_avg <- colMeans(island_clust1_TS)
island_clust2_avg <- colMeans(island_clust2_TS)

# island_x <- 1:length(island_clust1_avg)
# island_x <- years_kelp
island_x <- (1:length(island_clust1_avg) * 0.25) + 1983.75
# plot(main_x, mainland_clust1_avg, xlim=range(main_x), 
#      ylim=range(mainland_clust1_avg), xlab="Time", ylab="Kelp Biomass",
#      main="Mainland Average Time Series", pch=16, cex.lab=1.5, cex.main=1.8, cex.axis=1.8)

pdf(file=paste0(resloc, "island_avg_ts.pdf"), width=7, height=5)
plot(1, type="n", xlim=range(island_x), 
     ylim=range(island_clust2_avg), xlab="Year", ylab="Kelp biomass",
     # main="Time series of average kelp biomass: island locations", 
     cex.lab=1.7, cex.main=1.8, cex.axis=1.7,
     font.main=1)
lines(island_x[order(island_x)], island_clust1_avg[order(island_x)],
      xlim=range(island_x), ylim=range(island_clust1_avg), pch=16, col="#F8766D", lwd="2")
lines(island_x[order(island_x)], island_clust2_avg[order(island_x)],
      xlim=range(island_x), ylim=range(island_clust2_avg), pch=16, col="#00BFC4", lwd="2")
dev.off()

# for mainland data
mainland_clust1 <- which(main_cluster_numbers == 2)
mainland_clust2 <- which(main_cluster_numbers == 1)

mainland_clust1_TS <- main_cleanKelp[-c(mainland_clust2),]
mainland_clust2_TS <- main_cleanKelp[-c(mainland_clust1),]

mainland_clust1_avg <- colMeans(mainland_clust1_TS)
mainland_clust2_avg <- colMeans(mainland_clust2_TS)

# main_x <- 1:length(mainland_clust1_avg)
# main_x <- years_kelp
main_x <- (1:length(mainland_clust1_avg) * 0.25) + 1983.75
# plot(main_x, mainland_clust1_avg, xlim=range(main_x), 
#      ylim=range(mainland_clust1_avg), xlab="Time", ylab="Kelp Biomass",
#      main="Mainland Average Time Series", pch=16, cex.lab=1.5, cex.main=1.8, cex.axis=1.8)

pdf(file=paste0(resloc, "main_avg_ts.pdf"), width=7, height=5)
plot(1, type="n", xlim=range(main_x), 
     ylim=range(mainland_clust2_avg), xlab="Year", ylab="Kelp biomass",
     # main="Time series of average kelp biomass: mainland locations", 
     cex.lab=1.7, cex.main=1.8, cex.axis=1.7,
     font.main=1)
lines(main_x[order(main_x)], mainland_clust1_avg[order(main_x)],
      xlim=range(main_x), ylim=range(mainland_clust1_avg), pch=16, col="#F8766D", lwd="2")
lines(main_x[order(main_x)], mainland_clust2_avg[order(main_x)],
      xlim=range(main_x), ylim=range(mainland_clust2_avg), pch=16, col="#00BFC4", lwd="2")
dev.off()





