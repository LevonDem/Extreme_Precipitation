# Extreme_Precipitation
## Statistical modeling of extreme precipitation using TRMM data

This repository contains all of the R code used to generate the results in Demirdjian et al (2018), Statistical Modeling of Extreme Precipitation with TRMM Data. 

Before running any of the code, make sure you have downloaded the data (TRMM 3B42 daily data, spanning 1998/01/01 - 2014/02/04, in NetCDF format), available here: https://pmm.nasa.gov/data-access/downloads/trmm

### Step 1: Convert gridded data into individual time series (DataConversion.R).

The first step of this pipeline is to convert the gridded TRMM data, where each .nc data file is a grid containing the precipitation values for all grid points in a single day, into individual time series; each new file will correspond to the entire time series for a single grid point. This makes downstream analysis more efficient since we can then target a small subset of files corresponding to a specific region of interest.

This script takes a long time to run (around 1 day) and costs around 30GB of storage. 

To run this script, modify lines 31

``` setwd("~/Research/DataConversion") ```

and 36-37

``` 
main_dir   <- getwd()
setwd(paste(main_dir, "/netCDF_files", sep = ""))
```
to the directories containing the .nc files. Then simply source the script and let it run.

The file creates a directory for each longitude position (1-1440), each of which contains a separate file for each latitude position (1 - 400). 

### Step 2: Compute the 95th percentiles of precipitation levels at each grid point (MaxPrecipitation.R)

In the second step, we compute the 95th percentile of precipitation levels for each grid point. Of course, the script could easily be modified to compute a different statistic of interest, for example, the average precipitation at a given grid point. We will eventually use these computed percentiles when running our spatial clustering algorithm.

To run this script, simply modify the directories to reflect your setup:

Line 20:
``` main_dir <- getwd() ```

Line 33:
```   setwd(paste("~/DataConversion/", new.dir, sep = "")) ```

This script outputs a single text file, "Extreme_Prec_95.txt"

### Step 3: Spatial clustering of grid points (Clustering.R)

In the third step, we partition the 1440x400 grid points into non-overlapping clusters where extreme precipitation behavior can be considered to be reasonably homogeneous.

We use k-means clustering, performed in three iterations. We use the following covariates as input into the clustering algorithm:
1. Location (longitude and latitude).
2. 95th percentile of precipitation intensity. Alternatively, average precipitation can be used.
3. Topography (derived from 5' National Geophysical Data Center (NGDC) TerrainBase Global Digital Terrain Model (DTM), version 1.0, and binned into 0.25 resolution). This file is included in this repository.

The current version of the script runs the clustering algorithm using average precipitation as input. If you want to use the 90th/95th percentiles of precipitation instead (recommended), you can simply make the following change on line 46

``` 
avePrec <-  as.matrix(read.table("Extreme_Prec_95.txt", sep = " "))
```
to avoid renaming all of the variables in latter parts of the script. The output will be a .csv file, "cluster_results.csv", contanining the cluster number for each grid point in the TRMM map.

## Step 4: Compute average recurrence intervals (ARIs) over a target region (ARI_USA.R)

Finally, we fit a statistical model to each of the clusters in a specified region. In this script, we use the continental USA to illustrate the model fitting procedure. The fitted model is a non-stationary Poisson Point Process with month as a time covariate.

To use a different region, modify lines 102-106

```
long_small       <- 201:500
lat_small        <- 301:400
US_Clusters      <- clusterMap[long_small, lat_small]
US_Clusters_vec  <- unique(as.vector(US_Clusters))
numClusters      <- length(US_Clusters_vec)
```

To double check that this is the correct region, you can run the following plotting command:

```
my_plot(US_Clusters, long = longitude[long_small], lat = latitude[lat_small])
```

Note: the current version of the code computes confidence intervals for effective return levels. This adds alot of run time to the algorithm and can be skipped if we only want point estimates (for the USA, this version takes around 2.5 hours of runtime). 

In addition, this script returns the 2 and 25 year return levels. To change this, simply modify lines 268-269:
  
```  
rl.ci.NS1 <- try(ci(fit_MLE_NS, type = "return.level", qcov = v, return.period = 2))
rl.ci.NS2 <- try(ci(fit_MLE_NS, type = "return.level", qcov = v, return.period = 25))
```

The plotting commands near the end of the script should be modified to reflect the specific results you obtain (e.g. change the plot titles, limits, etc.)


