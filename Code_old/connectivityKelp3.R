# Miriam Wanner Connectivity and Kelp Biomass Data
# includes looking at points inside a certian band (more bidirectionally connected) - this part is taken out in 4th version of code
# April 12, 2021

library(spatstat)
library(geosphere)
library(wsyn)
library(ecodist)

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

# FIND SITES WITHIN A "BAND" variables:
# nSites = number of sites
# ptsProbMatFec = a matrix of four columns
#                 col 1: probMatFecundity as a vector
#                 col 2: transform of probMatFecundity as a vector
#                 col 3: the donor site
#                 col 4: the recipient site
# ptsInBand = same as ptsProbMatFec, except column 1, 2 for the ith row are NA when the ith row is not inside the band
# nEltsInBand = number of elements in the band
# ijInBand = column 1 is the donor site and column 2 is the recipient site. these are the points that are inside the band.

# find all the sites within a "band" of the probability matrix against itself
# this is when the dispersal from i to j is close to the dispersal from j to i
# matrix with two columns, each site's [i,j] dispersal and [j,i] dispersal
nSites <- nrow(probMatFecundity)
ptsProbMatFec <- matrix(c(as.vector(log10(probMatFecundity)), 
                          as.vector(log10(t(probMatFecundity))), 
                          rep(1:nSites, times = nSites), rep(1:nSites, each = nSites)), ncol = 4)
plot(ptsProbMatFec[,1], ptsProbMatFec[,2])
abline(a = 0, b = 1, col = "cyan")
abline(a = 3, b = 1, col = "red")
abline(a = -3, b = 1, col = "red")
# get the points inside certain lines, I chose y = x + 3 and y = x - 3, but I can change this
ptsInBand <- ptsProbMatFec
for(i in 1:nrow(ptsProbMatFec)){
  if(ptsProbMatFec[i, 3] == ptsProbMatFec[i, 4]){
    ptsInBand[i, 1] <- NA
    ptsInBand[i, 2] <- NA
  }
  else{
    ptx <- ptsProbMatFec[i, 1]
    pty <- ptsProbMatFec[i, 2]
    distLine1 <- abs(ptx - pty + 3) / sqrt(2)
    distLine2 <- abs(ptx - pty - 3) / sqrt(2)
    if(distLine1 > sqrt(18) || distLine2 > sqrt(18)){
      ptsInBand[i, 1] <- NA
      ptsInBand[i, 2] <- NA
    }
  }
}
plot(ptsInBand[,1], ptsInBand[,2])
abline(a = 0, b = 1, col = "cyan")
abline(a = 3, b = 1, col = "red")
abline(a = -3, b = 1, col = "red")

# make vector of i, j site locations
nEltsInBand <- length(ptsInBand[,1]) - length(which(is.na(ptsInBand[,1])))
ijInBand <- matrix(c(vector(,nEltsInBand), vector(,nEltsInBand)), ncol = 2)
k <- 1
for(i in 1:nrow(ptsInBand)){
  if(is.na(ptsInBand[i, 1]) == FALSE){
    ijInBand[k, 1] <- ptsInBand[i, 3]
    ijInBand[k, 2] <- ptsInBand[i, 4]
    k <- k + 1
  }
}

# SITES IN THE BAND SYNCHRONY AND PROBABILITY MATRICES variables:
# distInBand = the distances between the points in the band, all other points set to NA
# probInBand = probability*fecundity matrix, points not in band set to NA
# smInBand = sm from above, where points not in band set to NA
# symProbInBand = symProbMatF from above, where points not in band set to NA

# points where the sites i and j are most strongly connected in both directions
distInBand <- newDistMat
probInBand <- probMatFecundity
for(i in 1:nrow(ptsInBand)){
  if(is.na(ptsInBand[i, 1])){
    icoor <- ptsInBand[i, 3]
    jcoor <- ptsInBand[i, 4]
    distInBand[icoor, jcoor] <- NA
    probInBand[icoor, jcoor] <- NA
  }
}

plot(as.vector(distInBand), log10(as.vector(probInBand)))
plot(as.vector(distInBand[distInBand < 50]), log10(as.vector(probInBand[distInBand < 50])))

# synchrony where the sites i and j are most strongly connected in both directions
# (is it higher to places that are not as connected in both directions?)
smInBand <- sm
symProbInBand <- symProbMatF
for(i in 1:nrow(ptsInBand)){
  if(is.na(ptsInBand[i, 1])){
    icoor <- ptsInBand[i, 3]
    jcoor <- ptsInBand[i, 4]
    smInBand[icoor, jcoor] <- NA
    symProbInBand[icoor, jcoor] <- NA
  }
}
plot(log10(as.vector(smInBand)), log10(as.vector(symProbInBand)))

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


# USE MRM TO SEE CORRELATION BETWEEN SYNCHRONY AND SYMMETRIC PROBABILITY MATRIX variables:
MRM(sm ~ symProbMatF + newDistMat) # this gives me an error: "Matrix not square." but it is?



