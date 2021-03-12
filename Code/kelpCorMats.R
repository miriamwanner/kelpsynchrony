# March 5, 2021

library(wsyn)

# Setting up all the matrices
# Matrix for Biomass
file.choose()
kelpBioOriginal <- read.csv("/Users/miriam/Documents/Github/kelpsynchrony/Data/Kelp_Bio_2019_v3.csv")
# Matrix for Waves
file.choose()
wavesOriginal <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Hs_Max_2019.csv")
# Matrix for NO3 Concentration
file.choose()
NO3Original <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/NO3_Mean_2019.csv")


# Getting them all in paddle coordinates order
file.choose()
paddleCoordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Paddle_coords.csv")
order.paddle <- paddleCoordinates$Site_Number
# Biomass ordered by paddle coordinates
kelpBioPC <- kelpBioOriginal[match(order.paddle, kelpBioOriginal$Site),]
# Waves ordered by paddle coordinates
wavesPC <- wavesOriginal[match(order.paddle, wavesOriginal$Site),]
# NO3 ordered by paddle coordinates
NO3PC <- NO3Original[match(order.paddle, NO3Original$Site),]


# Getting only the coordinates from point conception to oxnard
# Using longitudes 120.471423 (point conception) to 119.16897 (oxnard)
# This is rows 263 to 339 (site numbers 275 to 339)
paddleCoorSB <- paddleCoordinates[263:339,]
kelpBioSB <- kelpBioPC[263:339,]
wavesSB <- wavesPC[263:339,]
NO3SB <- NO3PC[263:339,]


# Creating the correlation matrices
kelpBioCorNA <- cor(t(kelpBioSB), use="pairwise.complete.obs")
wavesBioCorNA <- cor(t(wavesSB), t(kelpBioSB), use="pairwise.complete.obs")
NO3BioCorNA <- cor(t(NO3SB), t(kelpBioSB), use="pairwise.complete.obs")
# Q: should I be transposing both the waves/NO3 matrix and the kelpBio matrix
# for the correlation matrix computation?


# Removing rows and columns that only contain NA
kelpBioCor <- kelpBioCorNA[rowSums(is.na(kelpBioCorNA)) != ncol(kelpBioCorNA), colSums(is.na(kelpBioCorNA)) != nrow(kelpBioCorNA)]
wavesBioCor <- wavesBioCorNA[rowSums(is.na(wavesBioCorNA)) != ncol(wavesBioCorNA), colSums(is.na(wavesBioCorNA)) != nrow(wavesBioCorNA)]
NO3BioCor <- NO3BioCorNA[rowSums(is.na(NO3BioCorNA)) != ncol(NO3BioCorNA), colSums(is.na(NO3BioCorNA)) != nrow(NO3BioCorNA)]


# Displaying the correlation matrices
image(kelpBioCor)
# places to note: 0.0 - 0.25 and 0.32 - 0.8 (indicate higher synchrony)
image(wavesBioCor)
# places to note (y axis): ~0.14, 0.25 - 0.3, 0.65 - 1
image(NO3BioCor)
# places to note (y axis): 0.25 - 0.3 (similar synchrony seen in waves here)
# Q: How are the graphs created from the waves/NO3 cor matrix?
# Why is there no symmetry for these?


# I wanted to try doing something with wavelets with the wsyn package, but 
# I do not think I understand wavelets enough yet. Is there a book I can 
# read to learn about them to be able to use the wsyn package?



