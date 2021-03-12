# March 10, 2021 Clustering and Wavelets

library(wsyn)

# Setting up all the matrices
# Matrix for Biomass
kelpBioOriginal <- read.csv("/Users/miriam/Documents/Github/kelpsynchrony/Data/Kelp_Bio_2019_v3.csv")

# Getting them all in paddle coordinates order
paddleCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
order.paddle <- paddleCoordinates$Site_Number
# Biomass ordered by paddle coordinates
kelpBioPC <- kelpBioOriginal[match(order.paddle, kelpBioOriginal$Site),]

# Getting only the coordinates from point conception to oxnard
# Using longitudes 120.471423 (point conception) to 119.16897 (oxnard)
# This is rows 263 to 339 (site numbers 275 to 339)
coords <- paddleCoordinates[263:339,]
kelpBioSB <- kelpBioPC[263:339,]

# Preparing the data for the 'clust' method
kelpBioSB <- subset(kelpBioSB, select=-Site) # remove the 'Site' column
# replace missing data with median of non missing values of row (in time series) - suggested in wsyn vignette
# is there a different/better method to get rid of the NA values?
medians <- apply(kelpBioSB, 1, median, na.rm=TRUE) # 1 for second arguments indicates rows
i <- 1 # is there an easier way than a for loop?
for(x in medians){ # for loop replaces each NA value in each row to be the median of that row
  if(is.na(x)){ # remove row when only contains NA in kelpBioSB data frame and coordinates
    kelpBioSB <- kelpBioSB[-i, ]
    coords <- coords[-i, ]
  }
  else{
    r <- kelpBioSB[i,]
    r[which(is.na(r))] <- x
    kelpBioSB[i,] <- r
    i <- i + 1
  }
}
times <- c(1:ncol(kelpBioSB)) # this is the times that will be passed into cleandat
matKelp <- data.matrix(kelpBioSB) # creates the matrix from the dataframe
cleanKelp <- cleandat(matKelp, times, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp)

# Synchrony matrix from the cleaned data - using pearson correlation
sm <- synmat(cleanKelp, 1:matcol, method="pearson")
fields::image.plot(1:ncol(sm),1:ncol(sm),sm,col=heat.colors(20)) # axes say "1:ncol(sm)" (change)

# Getting and plotting clusters
colnames(coords)[1] <- "lat" # rename columns for 'clust' method
colnames(coords)[2] <- "lon"
kelpClust <- clust(dat=cleanKelp,times=1:matcol,coords=coords,method="ReXWT")
get_clusters(kelpClust)
plotmap(kelpClust)
# creating wavelet mean fields for each module and plot
kelpClust <- addwmfs(kelpClust)
plotmag(get_wmfs(kelpClust)[[2]][[1]]) # I do not understand the indexing here but I just copied it from the vignette
plotmag(get_wmfs(kelpClust)[[2]][[2]]) 

# Wavelet transform for one location
kelpWT <- wt(t.series=cleanKelp[2,], times=times)
plotmag(kelpWT)
plotphase(kelpWT)


