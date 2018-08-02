##########################################################
## TRMM Precipitation data - Regional modeling in the USA
## Author: Levon Demirdjian
## Date : July 06, 2016
## Description: This script fits a separate distribution
## to each cluster in the continental USA
##
## Data: For each grid point in the TRMM daily data series 
## (1440 longitude and 400 latitude values), I created a 
## file containing the entire time series of daily 
## precipitation values, resulting in a total of 576,000 
## files. This script accesses these files.
##########################################################

library(doParallel) ## For parallel computation
library(extRemes)   ## For exreme value analysis
library(fields)     ## For plotting results
library(foreach)    ## For parallel computation


###############################
##  Set up plotting function ##
###############################

my_plot <- function(data, numColors = 100, long = longitude, lat = latitude, 
                    lower = min(data), upper = max(data), base.white = FALSE, 
                    ...){
  
  ## Define custom color scheme
  if(base.white == TRUE){
    colfunc <- colorRampPalette(c("white", "blue", "skyblue", "green", "yellow", "red"))
  }else{
    colfunc <- colorRampPalette(c("blue", "skyblue", "green", "yellow", "red"))
  }
  
  ## Set figure margins
  #par(mfrow=c(1,1))
  #par(mar = c(1, 1, 1, 1) + 0.1)
  
  ## Plot data 
  image.plot(long, lat, data, zlim = range(data),
             breaks = seq(lower, upper, length.out = numColors + 1),
             col = colfunc(numColors), xlab = "Longitude", ylab = "Latitude",
             xaxt = "n", yaxt = "n", ...)


  ## Overlay world map
  map("world", xlim = range(longitude), ylim = range(latitude), add = TRUE)
  
}


##############################################################
## Import precipitation, topography, and clustering results ##
##############################################################

## Import average precipitation data
avePrec    <- as.matrix(read.table("Average_Prec.txt", sep = " "))
topography <- as.matrix(read.table("Topography.txt", sep = " "))

## Import clustering results
clusterMap <- as.matrix(read.csv("cluster_results_07_28_2016.csv", header = TRUE, sep = ","))

## Define longitude and latitude
long_len   <- 1440
lat_len    <- 400
longitude  <- seq(-179.875, 179.875, length.out = long_len)
latitude   <- seq(-49.875, 49.875, length.out = lat_len)

#my_plot(clusterMap)

###################################
## Count number of days per year ##
###################################

main.dir   <- getwd()

## Set directory to wherever the raw netcdf files are 
setwd("~/netCDF_files")

## Count number of files per year
filenames  <- list.files(pattern = "*.nc")
years      <- 1998:2014
day_count  <- sapply(1:length(years), function(x) sum(grepl(years[x], filenames)))
setwd(main.dir)

## Organize dates according to the month
dates  <- read.csv("dates.csv", header = FALSE, sep = ",")

## Ignore the last 5 days since they are in 2014
month <- c()
for(i in 1:(length(filenames) - 5)){
  month[i] <- substr(as.character(dates[i,2]), 4, 6)
}


###################################
## Start model fitting procedure ##
###################################

## Only go over American Clusters
long_small       <- 201:500
lat_small        <- 301:400
US_Clusters      <- clusterMap[long_small, lat_small]
US_Clusters_vec  <- unique(as.vector(US_Clusters))
numClusters      <- length(US_Clusters_vec)

## Plot map for fun (optional)
#my_plot(US_Clusters, long = longitude[long_small], lat = latitude[lat_small])

## The cluster I used for diagnostic plots is 19339

## Setup parallel backend to use multiple processors
cl <- makeCluster(6)
registerDoParallel(cl)

## Create text file with run status
writeLines(c(""), "Algorithm_Status.txt")
start_time  <- Sys.time()
print(start_time)


## Store results into "results".
results <- foreach(iClust = US_Clusters_vec, .combine = cbind, .packages = c("extRemes")) %dopar% {
    
  ## Update status in the text file
  cur_index <- which(iClust == US_Clusters_vec)
  
  ## Monitor output while program runs
  sink("Algorithm_Status.txt", append=TRUE)
  cat(paste0(round(iClust / numClusters * 100), '% completed', '\n'), sep = "")
  sink()
  
  ## Find grid points corresponding to current cluster
  clustNumber <- iClust
  clust.ind   <- which(clusterMap == clustNumber, arr.ind = TRUE)
  
  ## Find average of average precipitation in this cluster
  clust.prec  <- mean(avePrec[clust.ind])
  
  ## Define months
  month_cov  <- rep(month, nrow(clust.ind)) 
  
  ## Import individual data
  main.dir   <- getwd()
  total.data <- list()
  for(i in 1:nrow(clust.ind)){
    
    ## Extract position of current grid point
    long.i  <- as.numeric(clust.ind[i,][1])
    lat.i   <- as.numeric(clust.ind[i,][2])
    
    ## Go into the relevant folder and extract the file corresponding to this grid point
    new.dir <- paste("Long", long.i, sep = "")
    setwd(paste("~/DataConversion/", new.dir, sep = ""))
    
    filename        <- paste("Long_", long.i, "_Lat_", lat.i, ".txt", sep = "")
    total.data[[i]] <- as.matrix(read.table(filename, sep = " "))[,1]
    total.data[[i]] <- total.data[[i]][1:5813]
    setwd(main.dir)
    
  }

  ## Pool all the data together and get rid of missing data
  is_valid    <- !is.na(unlist(total.data))
  pooled_data <- unlist(total.data)[is_valid]
  month_cov   <- month_cov[is_valid]
  
  #######################################################
  ##  Step 1: Find threshold and decluster pooled data ##
  #######################################################
  
  ## Find threshold for each month separately. There are 3 cases
  ## 1) If a region has very low average precipitation, set a lower
  ## bound and only consider precipitation levels above that bound.
  ## 2) If a region has only 1 grid point, use the 95th percentile
  ## to get some extra data values.
  ## 3) Else use the 99th percentile of the pooled series.
  
  monthly_thresh <- rep(0, nrow(clust.ind) * 5813)
  lower_bound    <- 0.1
  if(clust.prec <= lower_bound){
    thresh_quant    <- 0.95
    thresh          <- quantile(pooled_data[pooled_data >= lower_bound], probs = thresh_quant)
    monthly_thresh  <- rep(thresh, length(monthly_thresh))
  }else if((clust.prec > lower_bound) & (nrow(clust.ind) == 1)){
    thresh_quant    <- 0.95
    thresh          <- quantile(pooled_data, probs = thresh_quant)
    monthly_thresh  <- rep(thresh, length(monthly_thresh))
  }else{
    thresh_quant  <- 0.99
    thresh        <- quantile(pooled_data, probs = thresh_quant)
    
    ## Monthly thresholds
    thresh_ind_month <- list()
    for(i in 1:12){
      thresh_ind_month[[i]] <- quantile(pooled_data[month_cov == unique(month)[i]], probs = thresh_quant)
      monthly_thresh[month_cov == unique(month)[i]] <- thresh_ind_month[[i]]
    }
  }
  
  #############################
  ## SMOOTHING THE THRESHOLD ##
  #############################
  
  ## My choice of limits (e.g. 1904, 3000) is a bit ad hoc, but doesn't
  ## really make a big difference in the smoothing. This is to make sure
  ## that the function is smooth at the end points with no abrupt jumps
  spl_plot <- smooth.spline(1:1904, monthly_thresh[1097:3000], df = 90)
  
  ## Look at plot of smoothed threshold
  # plot(monthly_thresh, type = 'l', xlim = c(1, 365),
  #      xlab = "Day", ylab = "Precipitation",
  #      main = "Smoothed monthly threshold")
  # lines(spl_plot, col = "red", lwd = 3)
  
  monthly_spline_365 <- spl_plot$y
  monthly_spline <- c(monthly_spline_365[731:1856], 
                      monthly_spline_365[31:1856],
                      monthly_spline_365[31:1856],
                      monthly_spline_365[31:1065])
  
  final_thresh <- rep(monthly_spline, nrow(clust.ind))
  final_thresh <- final_thresh[is_valid]

  ax <- c()
  for(i in 1:nrow(clust.ind)){
    ax[i] <- sum(!is.na(total.data[[i]]))
  }
  ay <- c(0, cumsum(ax))
  
  ## Decluster individual series (using a window of size r = 5) and pool the results
  ind_series <- lapply(1:length(total.data), function(x) total.data[[x]][!is.na(total.data[[x]])])
  decl_data  <- unlist(lapply(1:length(total.data), function(x) decluster(ind_series[[x]], threshold = final_thresh[(1 + ay[x]):(ay[x+1])], r = 5)))
  n_size     <- sum(decl_data >= final_thresh)

  
  #######################################################
  ##  Step 2: Fit a model to declustered regional data ##
  #######################################################
  
  ## Fit declustered data for this region to a Point Process model.
  
  ## Define time covariates
  time_cov <- c(1:365, 1:365, 1:366, 1:365, 1:365, 1:365, 1:365, 
                1:366, 1:365, 1:365, 1:365, 1:365, 1:366, 
                1:365, 1:365, 1:335)
  time_cov  <- rep(time_cov, nrow(clust.ind))
  time_cov  <- time_cov[is_valid]

  ## Define non-stationary model with sinusoidal location and scale parameters
  fit_MLE_NS <- fevd(decl_data, threshold = final_thresh, type = "PP",
                     method = "MLE", use.phi = TRUE,
                     location.fun =~ sin(2 * pi * (time_cov)/365.25) +
                       cos(2 * pi * (time_cov )/365.25),
                     scale.fun =~  sin(2 * pi * (time_cov)/365.25) +
                       cos(2 * pi * (time_cov )/365.25))

  ## Prepare covariance matrix
  mu1  <- sin(2 * pi * (1:365) / 365.25)
  mu2  <- cos(2 * pi * (1:365) / 365.25)
  phi1 <- mu1
  phi2 <- mu2
  v    <- make.qcov(fit_MLE_NS, vals = list(mu1 = mu1, mu2 = mu2, phi1 = phi1, phi2 = phi2))
  
  ## Try computing confidence intervals for effective return levels; if error, return a constant
  ## This adds alot of run time to the algorithm and can be skipped if we only want point estimates.
  rl.ci.NS1 <- try(ci(fit_MLE_NS, type = "return.level", qcov = v, return.period = 2))
  rl.ci.NS2 <- try(ci(fit_MLE_NS, type = "return.level", qcov = v, return.period = 25))
  
  if("try-error" %in% class(rl.ci.NS1) | "try-error" %in% class(rl.ci.NS2)){
    rl.ci.NS1 <- matrix(1, nrow = length(mu1), ncol = 4)
    rl.ci.NS2 <- matrix(1, nrow = length(mu1), ncol = 4)
  }
  
  ## Save results
  result <- list(return_level_2 = rl.ci.NS1[1:365,],
                 return_level_25 = rl.ci.NS2[1:365,])
  result
  
}

## Close the clusters used for parallel computing
stopCluster(cl)
end_time  <- Sys.time()
print(end_time - start_time) ## 2.5 hours

## Save results
return_level_2  <- results[1,]
return_level_25 <- results[2,]


##################################
##  Plot 2 year ARI thresholds  ##
##################################

## Change day to whichever day of the year's ARI you want
day           <- 1
ARI_2         <- unlist(lapply(1:numClusters, function(x) return_level_2[[x]][day,2]))
ARI_2_map_Jan <- matrix(0, nrow = length(longitude), ncol = length(latitude))

## Takes a minute
for(i in 1:numClusters){
  clust.ind   <- which(clusterMap == US_Clusters_vec[i], arr.ind = TRUE)
  ARI_2_map_Jan[clust.ind] <- ARI_2[i]
}

## Save results for the United States
long_small     <- 200:500 
lat_small      <- 301:400 
ARI_2_map_sub_Jan  <- ARI_2_map_Jan[long_small, lat_small]

## Plot USA
my_plot(ARI_2_map_sub_Jan, 
        long = longitude[long_small], 
        lat = latitude[lat_small], 
        upper = 150, 
        main = paste('2 year ARI map on day', day, 'of the year'))
map("state", add = TRUE, lwd = 1)

## ARI_2_map_sub_Jan
## ARI_2_map_sub_July