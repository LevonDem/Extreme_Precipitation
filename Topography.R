########################################################
## TRMM Precipitation data - Compute Topography
## Author: Levon Demirdjian
## Date : July 22, 2016
## Description: This script computes the topography of 
## the world, zeroing values under sea level. We're 
## importing the file "topo_0.25deg_50S-50N.txt" here.
########################################################

## Switch directory to where the topographic data is
main_dir <- getwd()
setwd("~/TopographicData/data/")

## Import the data
top.data <- as.matrix(read.table("topo_0.25deg_50S-50N.txt"))

## Function to move all columns right; I'm converting
## longitude from (0, 360) to (-180, 180)
top_new <- matrix(0, nrow(top.data), ncol(top.data))
top_new[(((1:nrow(top.data))+719) %% 1440) + 1, ] <- top.data[1:nrow(top.data),]

## Set all negative topographic values to 0
top_new[top_new < 0] <- 0 

## Set directory back to what it was
setwd(main_dir)

## Write output as a .txt file
write.table(top_new, file = "Topography.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

## View topographic data (optional)
# longitude <- seq(-179.875, 179.875, length.out = 1440)
# latitude  <- seq(-49.875, 49.875, length.out = 400)
# 
# par(mar = c(2,3,2,2) + 0.1)
# image.plot(longitude, latitude, top_new, zlim = range(top_new),
#            breaks = seq(min(top_new), max(top_new), length.out = 100 + 1),
#            col = topo.colors(100, alpha = 1), xlab = "Longitude", 
#            ylab = "Latitude", main = "Topography")
# map("world", xlim = range(longitude), ylim = range(latitude), add = TRUE)
