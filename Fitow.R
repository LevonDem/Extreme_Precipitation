##########################################################
## TRMM Precipitation data - Modeling Typhoon Fitow
## Author: Levon Demirdjian
## Date : July 06, 2016
## Description: This script creates ARI maps for a
## a specific climate event, Typhoon Fitow.
##
## Data: For each grid point in the TRMM daily data series 
## (1440 longitude and 400 latitude values), I created a 
## file containing the entire time series of daily 
## precipitation values, resulting in a total of 576,000 
## files. This script accesses these files.
##########################################################

## 4 hours and 25 minutes to run
library(doParallel) ## For parallel computation
library(extRemes)   ## For exreme value analysis
library(fields)     ## For plotting results
library(foreach)    ## For parallel computation


###############################
##  Set up plotting function ##
###############################

my_plot <- function(data, numColors = 100, long = longitude, lat = latitude, 
                    lower = min(data), upper = max(data), base.white = FALSE, ...){
  
  ## Define custom color scheme
  if(base.white == TRUE){
    colfunc <- colorRampPalette(c("white", "blue", "skyblue", "green", "yellow", "red"))
  }else{
    colfunc <- colorRampPalette(c("blue", "skyblue", "green", "yellow", "red"))
  }
  
  ## Set figure margins
  #par(mfrow=c(1,1))
  par(mar = 3 * c(1, 1, 1, 1) + 0.1)
  
  ## Plot data 
  image.plot(long, lat, data, zlim = range(data),
             breaks = seq(lower, upper, length.out = numColors + 1),
             col = colfunc(numColors), xlab = "", ylab = "",
             xaxt = "n", yaxt = "n", ...)
  # xaxisLab <- paste(c(seq(180, 0, -30), seq(30, 180, 30)),
  #                   c(rep("W", 6), "", rep("E", 6)), 
  #                   sep = "")
  # 
  # yaxisLab <- paste(c(seq(50, 0, -10), seq(10, 50, 10)), 
  #                   c(rep("S", 5), "", rep("N", 5)), 
  #                   sep = "")
  # 
  # image(long, lat, data, zlim = range(data),
  #            breaks = seq(lower, upper, length.out = numColors + 1), 
  #            col = colfunc(numColors), xlab = "Longitude", 
  #            ylab = "Latitude", xaxt = "n", yaxt = "n",
  #            ...)
  # 
  # ## Add x-axis
  # axis(side = 1, at = seq(-180, 180, 30), labels = xaxisLab, 
  #      tck = 0)
  # 
  # ## Add y-axis
  # axis(side = 2, at = seq(-50, 50, 10), labels = yaxisLab, 
  #      las = 2, tck = 0)
  # 
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

## Typhoon Fitow info
FitowPrec  <- as.matrix(read.table("Fitow.txt", sep = " "))
long_small <- 1185:1220 
lat_small  <- 291:340 

## Look at clustering around the region where Fitow hit
clusterMapFitow <- clusterMap[long_small, lat_small]
#my_plot(clusterMapFitow, long = long_small, lat = lat_small)

###################################
## Count number of days per year ##
###################################

main.dir   <- getwd()
setwd("~/netCDF_files")

## Count number of files per year
filenames  <- list.files(pattern = "*.nc")
years      <- 1998:2014
day_count  <- sapply(1:length(years), function(x) sum(grepl(years[x], filenames)))
setwd(main.dir)

## Organize dates according to the month
dates  <- read.csv("dates.csv", header = FALSE, sep = ",")

## Ignore the last 5 days
month <- c()
for(i in 1:(length(filenames) - 5)){
  month[i] <- substr(as.character(dates[i,2]), 4, 6)
}


###################################
## Start model fitting procedure ##
###################################

FitowClusters <- unique(as.vector(clusterMap[long_small, lat_small]))
numClusters   <- length(FitowClusters)

## Setup parallel backend to use multiple processors
cl <- makeCluster(6)
registerDoParallel(cl)

## Create text file with run status
writeLines(c(""), "Algorithm_Status.txt")
start_time  <- Sys.time()
print(start_time)


## Store results into "results".
results <- foreach(iClust = FitowClusters, .combine = cbind, .packages = c("extRemes")) %dopar% {

  ## Update status in the text file
  cur_index <- which(iClust == FitowClusters)
  
  sink("Algorithm_Status.txt", append=TRUE)
  cat(paste0(round(cur_index * 100 / 270), '% completed', '\n'), sep = "")
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
  
  ## SMOOTHING THE THRESHOLD
  ## Only smooth days 1:1100, then repeat smoothed series
  spl_plot <- smooth.spline(1:1904, monthly_thresh[1097:3000], df = 90)
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
  
  ## Label these more intuitively
  ax <- c()
  for(i in 1:nrow(clust.ind)){
    ax[i] <- sum(!is.na(total.data[[i]]))
  }
  ay <- c(0, cumsum(ax))
  
  ## Decluster individual series and pool the results
  ind_series <- lapply(1:length(total.data), function(x) total.data[[x]][!is.na(total.data[[x]])])
  decl_data  <- unlist(lapply(1:length(total.data), function(x) decluster(ind_series[[x]], threshold = final_thresh[(1 + ay[x]):(ay[x+1])], r = 5)))
  n_size     <- sum(decl_data >= final_thresh)
  
  #######################################################
  ##  Step 2: Fit a model to declustered regional data ##
  #######################################################
  
  ## Fit declustered data for this region to a GP distribution
  ## We can use L-moments or Bayesian as well, but the sample 
  ## size appears to be sufficiently large that the results are 
  ## very similar to one another.
  
  time_cov <- c(1:365, 1:365, 1:366, 1:365, 1:365, 1:365, 1:365, 
                1:366, 1:365, 1:365, 1:365, 1:365, 1:366, 
                1:365, 1:365, 1:335)
  time_cov  <- rep(time_cov, nrow(clust.ind))
  time_cov  <- time_cov[is_valid]
  
  ## Define non-stationary model
  fit_MLE_NS <- fevd(decl_data, threshold = final_thresh, type = "PP",
                     method = "MLE", use.phi = TRUE,
                     location.fun =~ sin(2 * pi * (time_cov)/365.25) +
                       cos(2 * pi * (time_cov )/365.25),
                     scale.fun =~  sin(2 * pi * (time_cov)/365.25) +
                       cos(2 * pi * (time_cov )/365.25))


  ## Compute return levels
  year_list <- c(1.01, 2:200)

  ## Prepare covariance matrix
  FitowDay <- 279
  mu1      <- sin(2 * pi * (FitowDay) / 365.25)
  mu2      <- cos(2 * pi * (FitowDay) / 365.25)
  phi1     <- mu1
  phi2     <- mu2
  v        <- make.qcov(fit_MLE_NS, vals = list(mu1 = mu1, mu2 = mu2, phi1 = phi1, phi2 = phi2))

  rl.ci.NS <- list()
  for(i in 1:100){
    rl.ci.NS[[i]] <- ci(fit_MLE_NS, type = "return.level", qcov = v, return.period = year_list[i])
  }
  
  ## Save results
  result <- list(return_level = rl.ci.NS)
  result
  
}

## Close the clusters used for parallel computing
stopCluster(cl)
end_time  <- Sys.time()
print(end_time - start_time) ## 2.5 hours

## Save results
return_levels <- results[1,]
rm(results)

############################
## Store results onto map ##
############################

ARI_Fitow <- matrix(0, nrow = nrow(FitowPrec), ncol = ncol(FitowPrec))
for(i in 1:nrow(clusterMapFitow)){
  for(j in 1:ncol(clusterMapFitow)){
    
    ## Find which cluster this point belongs to
    clust.ind <- which(FitowClusters == clusterMapFitow[i, j])

    ## Get return values for that cluster
    rl_ind <- unlist(lapply(1:100, function(x) return_levels[[clust.ind]][[x]][1,2]))
    
    ## Find which rl year is closest to observed precipitation
    best_rl <- which.min(abs(FitowPrec[i, j] - rl_ind))
    
    ## Save year on map
    ARI_Fitow[i, j] <- best_rl
    
    # if(best_rl == 100)
    #   cat(clusterMapFitow[i, j], "\n")
  }
}
  

## Plot precipitation map
FitowPrecWhite <- FitowPrec
FitowPrecWhite[which(FitowPrecWhite <= 100, arr.ind = TRUE)] <- 10000
my_plot(FitowPrecWhite, long = longitude[long_small], lat = latitude[lat_small],
        numColors = 300, lower = 50, upper = 350,
        main = paste("Precipitation (6 Oct 2013)"))
map.axes(cex.axis=1)
mtext("[mm]", at = c(127.15,30), cex = 0.85)

## Plot ARI map
FitowARIWhite <- ARI_Fitow
FitowARIWhite[which(FitowPrecWhite >= 1000, arr.ind = TRUE)] <- 1000
my_plot(FitowARIWhite, long = longitude[long_small], lat = latitude[lat_small],
        numColors = 25, lower = 0, upper = 100,
        main = paste("ARI (6 Oct 2013)"))
map.axes(cex.axis=1)
mtext("[year]", at = c(127.25,30), cex = 0.85)
mtext("+", at = c(127.60,19), cex = 0.85, line = -1.25)
