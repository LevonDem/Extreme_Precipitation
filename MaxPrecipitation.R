#############################################################
## TRMM Precipitation data - Extreme precipitation
## Author: Levon Demirdjian
## Date : July 26, 2016
## Description: This script computes the 95th (or 90th) 
## percentiles of precipitation levels at each grid 
## point for the TRMM. 
##
## Data structure: The data this file accesses is organized
## in the following form:
## There are 1440 folders corresponding to each longitude
## value (named Long1, Long2, etc.),  each of which has 400
## files corresponding to each latitude value. The names of
## these files are, for example, Long_1_Lat_1.txt, 
## Long_1_Lat_2.txt, etc. Each file contains the time series
## for the entire TRMM data record at that grid point.
#############################################################

## Save current directory
main_dir <- getwd()

## Define the resolution of the data images
numLong     <- 1440
numLat      <- 400
percentiles <- matrix(0, nrow = numLong, ncol = numLat)

## Import all of the data
for(i in 1:numLong){
  cat("We're on iteration ", i, " out of ", numLong, ".\n", sep = "")
  
  ## Change directory to current longitude's folder
  new.dir <- paste("Long", i, sep = "")
  setwd(paste("~/DataConversion/", new.dir, sep = ""))
  for(j in 1:numLat){
    
    ## Import current latitude's data series
    filename  <- paste("Long_", i, "_Lat_", j, ".txt", sep = "")
    ind.data  <- as.matrix(read.table(filename, sep = " "))[,1]
    
    ## Eliminate missing values
    ind.data  <- ind.data[!is.na(ind.data)]
    
    ## Save percentiles (I'm using the 95 percentile here)
    percentiles[i, j] <- as.numeric(quantile(ind.data, probs = 0.95))
  }
  
  ## Move directory up one level
  setwd(main_dir)
}


## Write output as a .txt file
write.table(percentiles, file = "Extreme_Prec_90.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

## View precipitation percentile data (optional)
# longitude <- seq(-179.875, 179.875, length.out = 1440)
# latitude  <- seq(-49.875, 49.875, length.out = 400)
# 
# par(mar = c(2,3,2,2) + 0.1)
# image.plot(longitude, latitude, percentiles, zlim = range(percentiles),
#            breaks = seq(min(percentiles), max(percentiles), length.out = 100 + 1),
#            col = topo.colors(100, alpha = 1), xlab = "Longitude", 
#            ylab = "Latitude", main = "Precipitation percentiles")
# map("world", xlim = range(longitude), ylim = range(latitude), add = TRUE)
