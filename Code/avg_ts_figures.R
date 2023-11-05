# creating the average time series figures for the clusters

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

load(file=paste0(resloc, "clustering_figures.RData"))

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

all_x <- (1:length(all_clust1_avg) * 0.25) + 1983.75

pdf(file=paste0(resloc, "all_avg_ts.pdf"), width=7, height=5)
par(mar=c(5,5,1,1))
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

island_x <- (1:length(island_clust1_avg) * 0.25) + 1983.75

pdf(file=paste0(resloc, "island_avg_ts.pdf"), width=7, height=5)
par(mar=c(5,5,1,1))
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

main_x <- (1:length(mainland_clust1_avg) * 0.25) + 1983.75

pdf(file=paste0(resloc, "main_avg_ts.pdf"), width=7, height=5)
par(mar=c(5,5,1,1))
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



