# creating the map clustering figures

library("here")
here::i_am("kelpsynchrony/Code/main.R")
setwd(here("kelpsynchrony/Code"))
rm(list=ls())
# location for storing the results
resloc <- "../Results/"
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
# ORGANIZE ROMS SPATIAL DATA

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

# save(roms.pairs.list, file=here::here("roms.pairs.list.rda"))

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


# clustering map

times <- c(1:ncol(kelpDataROMSSites))

colnames(newLonLatROMS)[2] <- "lat"
colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
coords <- as.data.frame(newLonLatROMS)

cleanKelp <- cleandat(data.matrix(kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp)
kelpClust <- clust(dat=cleanKelp,times=1:matcol,coords=coords,method="pearson")
cluster_numbers <- get_clusters(kelpClust)[[2]]

roms.pairs.list.subset <- data.frame(kelpClust$coords$lon, kelpClust$coords$lat, kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


coastline.data.dir        <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_f.b")
coastline.border.data.dir <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_f.b")
inset.data.dir            <- paste0(datloc, "gshhg-bin-2.3.7/gshhs_l.b")
inset.border.data.dir     <- paste0(datloc, "gshhg-bin-2.3.7/wdb_borders_l.b")

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

# Set parameters for maps
lat.0 <- 32.45
lat.1 <- 34.65

long.0 <- -120.65
long.1 <- -117


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




# clustering with all data

roms.pairs.list.subset <- data.frame(kelpClust$coords$lon, kelpClust$coords$lat, kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "all_clustering_map.pdf"), width=7, height=5)
all.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
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
all.map
dev.off()


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

roms.pairs.list.subset <- data.frame(island_kelpClust$coords$lon, island_kelpClust$coords$lat, island_kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "island_clustering_map.pdf"), width=7, height=5)
island.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
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
island.map
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


roms.pairs.list.subset <- data.frame(main_kelpClust$coords$lon, main_kelpClust$coords$lat, main_kelpClust$clusters[[2]])
names(roms.pairs.list.subset)[1] <- "lon"
names(roms.pairs.list.subset)[2] <- "lat"
names(roms.pairs.list.subset)[3] <- "Cluster"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 1] <- "Cluster 1"
roms.pairs.list.subset$Cluster[roms.pairs.list.subset$Cluster == 2] <- "Cluster 2"


pdf(file=paste0(resloc, "mainland_clustering_map.pdf"), width=7, height=5)
mainland.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=lon, y=lat, fill=Cluster, color=Cluster), size=2.2, shape=21) + # , alpha=0.5) +
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
mainland.map
dev.off()


save.image(file = paste0(resloc, "clustering_figures.RData"))

