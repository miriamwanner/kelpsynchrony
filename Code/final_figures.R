# script for final figures

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/ROMSKelpData/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"
load(file=paste0(resloc, "intermediate_figures.RData"))


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

library(dplyr)
library(tidyr)
library(reshape2)
library(sp)
library(rgdal)
library(maptools)
library(ggsn)
library(patchwork)


# ========================== 8 models with MRM ==========================

# USE MRM TO SEE CORRELATION BETWEEN SYNCHRONY AND SYMMETRIC PROBABILITY MATRIX:


MRM(as.dist(sm) ~ as.dist(logSym) + as.dist(newDistMat) + as.dist(smNO3))
MRM(as.dist(sm) ~ as.dist(logSym))

MRM(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))

MRM(as.dist(sm) ~ as.dist(logSym) + as.dist(newDistMat))
MRM(as.dist(sm) ~ as.dist(logSym) + as.dist(smNO3))

MRM(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(smNO3))


MRM(as.dist(sm) ~ as.dist(newDistMat) + as.dist(logSym))
MRM(as.dist(sm) ~ as.dist(logSym))
MRM(as.dist(sm) ~ as.dist(newDistMat))




# logit transform
diag(sm) <- NA
logit_sm <- logit(sm, min=-1, max=1)
diag(sm) <- 1

# vector_logit_sm = logit_sm[upper.tri(logit_sm, diag=FALSE)]
# logit_max = max(vector_logit_sm)
# diag(logit_sm) <- logit_max

# linear - transport: symProbMat, connectivity: symProbMatF
# log - transport: logSymProbMat, connectivity: logSym
# logit: logit_sm, not logit: sm

MRM(as.dist(sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))

source("altered_mrm_function.R")
# anytime i call MRM, also call my_mrm and make sure the results are the same
my_mrm(as.dist(sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))





# ========================== Spline Correlogram ==========================


# spline correlogram


spline_corr = Sncf(newLonLatROMS[,1], newLonLatROMS[,2], kelpDataROMSSites, latlon = TRUE)
synch_vec <- sm[upper.tri(sm, diag=FALSE)]
dist_vec <- newDistMat[upper.tri(newDistMat, diag=FALSE)]
smoothingSpline = smooth.spline(spline_corr$real$predicted$x, spline_corr$real$predicted$y, spar=0.35)
plot(spline_corr, ylim=c(-0.25, 1), cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
points(dist_vec, synch_vec, pch=20, cex=0.2, col="lightgrey")
lines(smoothingSpline)
# points(spline_corr$real$predicted$x, spline_corr$real$predicted$y, pch=20, cex=0.2)
# spline_corr$boot$boot.summary$predicted has different rows

# ========================== Clustering Map ==========================

times <- c(1:ncol(kelpDataROMSSites))

colnames(newLonLatROMS)[1] <- "lon" # rename columns for 'clust' method
colnames(newLonLatROMS)[2] <- "lat"
coords <- as.data.frame(newLonLatROMS)

cleanKelp <- cleandat(data.matrix(kelpDataROMSSites), times, clev=5)$cdat # cleans the data - max cleaning level
matcol <- ncol(cleanKelp)
kelpClust <- clust(dat=cleanKelp,times=1:matcol,coords=coords,method="pearson")
cluster_numbers <- get_clusters(kelpClust)[[3]]
plotmap(kelpClust)

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
coast.polys.sp <- PolySet2SpatialPolygons(coast.polys, close_polys=FALSE) # Convert to spatial polygons
polys <- fortify(coast.polys)

# Coastline borders
coast.border <- importGSHHS(gshhsDB=coastline.border.data.dir,
                            xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.border$X <- coast.border$X - 360  # Get longitude back into units of degrees east
coast.border.sp <- PolySet2SpatialPolygons(coast.border, close_polys=FALSE) # Convert to spatial polygons
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
                 st.bottom = FALSE, height=0.02, st.dist = 0.035, st.size = 6,
                 st.color = 'black', box.fill = 'black', box.color = 'transparent') +
  scale_x_continuous(breaks=rev(seq(floor(long.0), ceiling(long.1), by = 1)), 
                     labels= -1 * rev(seq(floor(long.0), ceiling(long.1), by = 1))) +
  scale_y_continuous(breaks= seq(floor(lat.0), ceiling(lat.1), by=1),
                     labels= seq(floor(lat.0), ceiling(lat.1), by = 1)) +
  xlab(expression(paste(Longitude," (",degree,W,")", sep=""))) +
  ylab(expression(paste(Latitude," (",degree,N,")", sep=""))) +
  annotate(geom = 'text', x = -119.75, y = 32.9, label = 'Pacific Ocean', size = 6) +
  annotate(geom = 'text', x = -119.25, y = 34.6, label = 'California, USA', size = 6)+
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


# getting the different colors of clusters
roms.pairs.list.subset$Region <- rep(NA, 43)
roms.pairs.list.subset$Region[which(cluster_numbers == 1)] <- c(rep("Cluster 1", length(which(cluster_numbers == 1))))
roms.pairs.list.subset$Region[which(cluster_numbers == 2)] <- c(rep("Cluster 2", length(which(cluster_numbers == 2))))
roms.pairs.list.subset$Region[which(cluster_numbers == 3)] <- c(rep("Cluster 3", length(which(cluster_numbers == 3))))
# roms.pairs.list.subset$i.roms.lat <- newLonLatROMS[,2]
# roms.pairs.list.subset$i.roms.long <- newLonLatROMS[,1]


roms.map <- base.map + 
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat, fill=Region), size=9, shape=21, alpha=0.5) +
  geom_point(data=roms.pairs.list.subset, aes(x=i.roms.long, y=i.roms.lat), size=1, color="black") +
  geom_point(data=kelpCoordinates, aes(x=Lon, y=Lat), size=0.5, color="grey") +
  theme_bw() + theme(# Set custom formatting
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", size=1.5),
    panel.background = element_blank()) + 
  theme(legend.key = element_blank(), legend.background=element_blank()) +
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



# ========================== Colorful predictor matrices ==========================


# making mulitpanel plots
par(mfrow=c(3,2), mai = c(0.3, 0.3, 0.3, 0.3))
# layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
# colorful matrix of sm
fields::image.plot(1:ncol(sm),1:ncol(sm),sm,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of smWaves
fields::image.plot(1:ncol(smWaves),1:ncol(smWaves),smWaves,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5) # axes say "1:ncol(sm)" (change)
# colorful matrix of smNO3
fields::image.plot(1:ncol(smNO3),1:ncol(smNO3),smNO3,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5) # axes say "1:ncol(sm)" (change)
# colorful matrix of newDistMat
fields::image.plot(1:ncol(newDistMat),1:ncol(newDistMat),newDistMat,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of transport - log logSymProbMat
diag(logSymProbMat) <- 0
fields::image.plot(1:ncol(logSymProbMat),1:ncol(logSymProbMat),logSymProbMat,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
# colorful matrix of connectivity - log logSym
fields::image.plot(1:ncol(logSym),1:ncol(logSym),logSym,col=heat.colors(20),xlab="Sites",ylab="Sites",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)


save.image(file = paste0(resloc, "final_figures.RData"))



