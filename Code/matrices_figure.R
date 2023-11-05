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



# colorful matrices


pdf(file=paste0(resloc, "matrices_figure.pdf"), width=7, height=10)
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
