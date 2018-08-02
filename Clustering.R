########################################################
## TRMM Precipitation data - Multilevel clustering
## Author: Levon Demirdjian
## Date : July 05, 2016
## Description: This script carries out multiple 
## stages of k-means clustering, starting from a coarse 
## clustering of the map, then further clustering each 
## of the resulting regions.
########################################################


###############################
##  Load necessary libraries ##
###############################

library(fields)   ## For plotting results

###############################
##  Set up plotting function ##
###############################

my_plot <- function(data, numColors, ...){
  
  ## Define custom color scheme
  colfunc <- colorRampPalette(c("blue", "skyblue", "green", "yellow", "red"))
  
  ## Set figure margins
  par(mar = c(1, 1, 1, 1) + 0.1)
  
  ## Plot data 
  image.plot(longitude, latitude, data, zlim = range(data),
             breaks = seq(min(data), max(data), length.out = numColors + 1), 
             col = colfunc(numColors), xlab = "Longitude", ylab = "Latitude", ...)
  
  ## Overlay world map
  map("world", xlim = range(longitude), ylim = range(latitude), add = TRUE)

}


#############################################################
##  Step 1: Import data to be used in clustering algorithm ##
#############################################################

## Import average precipitation and 90th percentiles of precipitation
avePrec     <- as.matrix(read.table("Average_Prec.txt", sep = " "))
extremePrec <- as.matrix(read.table("Extreme_Prec_90.txt", sep = " "))

## Import topography data
topography <- as.matrix(read.table("Topography.txt", sep = " "))


####################################
##  Step 2: High-level clustering ##
####################################

## Define longitude and latitude grids
numLong   <- 1440
numLat    <- 400
longitude <- seq(-179.875, 179.875, length.out = numLong)
latitude  <- seq(-49.875, 49.875, length.out = numLat)

## Vectorize the data for clustering
loc_scale   <- 1 ## Weights to give longitude and latitude (can change this)
avePrecVec  <- as.vector(avePrec)
topogVec    <- as.vector(topography)
long_vec    <- rep(longitude,  length(latitude))       
lat_vec     <- rep(latitude, each = length(longitude))

## Standardize the data for clustering
## Note: We need to wrap aroud the map, hence why I'm
## copying things over twice.
avePrecVec  <- scale(c(avePrecVec, avePrecVec))
topogVec    <- scale(c(topogVec, topogVec))
long_vec    <- loc_scale * scale(1:length(avePrecVec)) ## For padding the map
lat_vec     <- loc_scale * scale(c(lat_vec, lat_vec))

## Combine all the data that is going to be used for clustering
clust_data  <- data.frame(avePrecVec, topogVec, long_vec, lat_vec)

## Perform k-means clustering
numClusters  <- 60
start.time   <- Sys.time()

## Start clustering: this takes around 8 minutes for 60 clusters, 10 starts
clust_res    <- kmeans(clust_data, centers = numClusters, algorithm ="MacQueen", iter.max = 1000, nstart = 10) 

## Since there is wrap around, we keep only the first half of the results
clusterMap_t <- matrix(clust_res$cluster[1:(length(avePrecVec)/2)], length(longitude), length(latitude))
end.time     <- Sys.time()
print(end.time - start.time) 

## Reorganize cluster values to be consecutive starting from 1
numClusters <- length(unique(as.vector(clusterMap_t)))
clusterMap  <- matrix(0, nrow(clusterMap_t), ncol(clusterMap_t))
for(i in 1:numClusters){
  clust.ind             <- which(clusterMap_t == sort(unique(as.vector(clusterMap_t)))[i])
  clusterMap[clust.ind] <- i
}

## Plot high level clustering results (optional)
# my_plot(clusterMap, numColors = numClusters, main = "High level clusters")


###################################
##  Step 3: Mid level clustering ##
###################################

numClusters2  <- 30
clusterMap2   <- matrix(0, length(longitude), length(latitude))
start.time    <- Sys.time()
for(iCluster in 1:numClusters){
  
  loc_scale2   <- 1 ## Weights for longitude and latitude
  clust_index  <- which(clusterMap == iCluster, arr.ind = TRUE)
  
  ## Vectorize and scale the data
  avePrecVec2   <- scale(as.vector(avePrec[clust_index]))
  topogVec2     <- scale(as.vector(topography[clust_index]))
  long_vec2     <- loc_scale2 * scale(longitude[clust_index[,1]])       
  lat_vec2      <- loc_scale2 * scale(latitude[clust_index[,2]])
  
  ## If there is no variation in topography, do not use it as a 
  ## covariate for clustering
  if(sum(!is.na(topogVec2)) == 0){
    clust_data2  <- data.frame(avePrecVec2, long_vec2, lat_vec2)
  }else{
    clust_data2  <- data.frame(avePrecVec2, topogVec2, long_vec2, lat_vec2)
  }
  
  ## Perform k-means clustering if there are enough grid points in cluster
  if(nrow(clust_index) > numClusters2){
    clust_res2 <- kmeans(clust_data2, centers = numClusters2, algorithm ="MacQueen", iter.max = 1000, nstart = 10) 
    clusterMap2[clust_index] <- max(clusterMap2) + clust_res2$cluster
  }else{
    clusterMap2[clust_index] <- max(clusterMap2) + 1
  }

}
end.time     <- Sys.time()
print(end.time - start.time) ## 50 seconds for 30 clusters

## Plot mid level clustering results (optional)
# my_plot(clusterMap2, numColors = length(unique(as.vector(clusterMap2))), main = "Mid level clusters")


####################################
##  Step 4: Fine level clustering ##
####################################

numClusters3 <- 30
clusterMap3  <- matrix(0, length(longitude), length(latitude))
start.time   <- Sys.time()
for(iCluster in 1:(numClusters * numClusters2)){
  
  loc_scale3   <- 1 ## Weights for longitude and latitude
  clust_index  <- which(clusterMap2 == iCluster, arr.ind = TRUE)
  
  ## Vectorize and scale the data
  topogVec3     <- scale(as.vector(topography[clust_index]))
  avePrecVec3   <- scale(as.vector(avePrec[clust_index]))
  long_vec3     <- loc_scale3 * scale(longitude[clust_index[,1]])         
  lat_vec3      <- loc_scale3 * scale(latitude[clust_index[,2]])
  
  ## If there is no variation in topography, do not use it as a 
  ## covariate for clustering
  if(sum(!is.na(topogVec3)) == 0){
    clust_data3  <- data.frame(avePrecVec3, long_vec3, lat_vec3)
  }else{
    clust_data3  <- data.frame(avePrecVec3, topogVec3, long_vec3, lat_vec3)
  }
  
  ## Perform k-means clustering if there are enough grid points in cluster
  if(nrow(clust_index) > numClusters3){
    clust_res3    <- kmeans(clust_data3, centers = numClusters3, algorithm ="MacQueen", iter.max = 1000, nstart = 10) 
    clusterMap3[clust_index] <- max(clusterMap3) + clust_res3$cluster
  }else{
    clusterMap3[clust_index] <- max(clusterMap3) + 1
  }
  
}
end.time <- Sys.time()
print(end.time - start.time) ## 24 seconds

## Plot fine level clustering results (optional)
numClustersTotal <- length(unique(as.vector(clusterMap3)))
my_plot(clusterMap3, numColors = numClustersTotal, main = "Clustering results")

## Write clustering output to a .csv file
# write.csv(clusterMap3, file = "cluster_results_07_28_2016.csv", quote = FALSE, sep = " ",
#           row.names = FALSE, col.names = FALSE)

