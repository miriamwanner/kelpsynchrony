# February 19, 2021

library("ggcorrplot")
library("tidyverse")
library("ggmap")
library("ggplot2")
library("sf")
library("mapview")
library("maps")
library("mapdata")

file.choose()
KelpBio <- read.csv("/Users/miriam/Documents/Github/kelpsynchrony/Data/Kelp_Bio_2019_v3.csv")
View(KelpBio)

# Correlation Matrices
# kelpCor uses "pairwise" and kelpCor2 uses "complete.obs" to address NaN values
kelpCor <- cor(KelpBio, use="pairwise") # gives error that standard deviation is zero
kelpCor2 <- cor(KelpBio, use="complete.obs") # gives error that standard deviation is zero
View(kelpCor)
View(kelpCor2)
p <- ggcorrplot(kelpCor)
p + theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
  ggtitle("Correlation with pairwise")
p2 <- ggcorrplot(kelpCor2)
p2 + theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
  ggtitle("Correlation with complete.obs")

# Plot Data Locations on a Map (I tried this but it didn't work because I need an API to use google maps)
file.choose()
coordinates <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Sites.csv")
View(coordinates)
coordinates_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = 4326)
mapview(coordinates_sf)
#map <- map_data("california")
#ggplot() + 
#  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = NA, color = "red") + 
#  coord_fixed(1.3) + geom_point(data = coordinates, aes, x = Lon, y = Lat)
#map <- get_googlemap(center = c(-119.7, 34.4), zoom = 6,
#                     color = "bw",
#                     style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
#ggmap(map) + geom_point(data = coordinates, aes, x = Lon, y = Lat)

# How wave data correlates to biomass for a specific location
file.choose()
waves <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/Hs_Max_2019.csv")
View(waves)
wave1 <- as.vector(waves[4,])
biomass1 <- as.vector(KelpBio[4,])
waveBioData <- data.frame(t(wave1), t(biomass1))
View(waveBioData)
waveAndBio <- cor(waveBioData, use="complete.obs")
View(waveAndBio)

# How wave height data correlates to biomass for more locations
allWavesAndBio <- cor(waves, KelpBio, use="complete.obs")
View(allWavesAndBio)
w <- ggcorrplot(allWavesAndBio)
w + theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
  ggtitle("Correlation of wave heights and biomass")

# How NO3 concentration data correlates to biomass for all locations
file.choose()
NO3 <- read.csv("/Users/miriam/Documents/GitHub/kelpsynchrony/Data/NO3_Mean_2019.csv")
no3AndBio <- cor(NO3, KelpBio, use="complete.obs")
View(no3AndBio)
n <- ggcorrplot(no3AndBio)
n + theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
  ggtitle("Correlation of NO3 and biomass")


# other things I could do
# - synchrony (I do not know how I would calculate this)
# - using wavelets (like in the papers I read)
# - mapping synchrony (like in geography of synchrony paper)

