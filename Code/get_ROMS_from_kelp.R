# Script for getting the ROMS sites from the kelp data

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/ROMSKelpData/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"


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
library(mms)


# PREPARING KELP DATA:
# kelpBioOriginal = original kelp data from CSV file
# paddleCoordinates = original paddle coordinates from CSV file
# order.paddle = vector of the sites in order of paddle coordinates
# kelpBioPC = kelp data with sites ordered by paddle coordinate
# whichSitesNA = a vector listing which sites were removed (because they were all NA)

# Matrix for Kelp Biomass
kelpBioOriginal <- read.csv(file=paste0(datloc, "Kelp_Bio_2019_v3.csv"))
paddleCoordinates <- read.csv(file=paste0(datloc, "Paddle_coords.csv"))
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

OceanDistanceMatrices <- readMat(con=paste0(datloc, "OceanDistanceMatrices.mat"))
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
kelpCoordinates <- read.csv(file=paste0(datloc, "Paddle_coords.csv"))
kelpCoordinates <- kelpCoordinates[-whichSitesNA, ] # remove the sites that were removed earlier (because the entire time series was NA)
romsCoordinates <- as.data.frame(t(readMat(con=paste0(datloc, "site_centers.mat")))) # this is the center of the coordinates
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

roms.pairs.list.subset <- subset(roms.pairs.list.mainland.SN[[1]], !is.na(ROMSSitesVec))

save.image(file = paste0(resloc, "get_ROMS.RData"))


# save(order.paddle, file = paste0(resloc, "order.paddle.RData"))
# save(whichSitesNA, file = paste0(resloc, "whichSitesNA.RData"))
# save(romsCoordinates, file = paste0(resloc, "romsCoordinates.RData"))
# save(vecLocations, file = paste0(resloc, "vecLocations.RData"))
# save(newSites, file = paste0(resloc, "newSites.RData"))
# save(kelpDataROMSSites, file = paste0(resloc, "kelpDataROMSSites.RData"))
# save(oceanAvgKelpSites, file = paste0(resloc, "oceanAvgKelpSites.RData"))
