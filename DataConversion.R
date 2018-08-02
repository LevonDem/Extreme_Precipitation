############################################################
## TRMM Precipitation data - Data conversion
## Author: Levon Demirdjian
## Date : July 08, 2016
## Description: This script converts the gridded TRMM data,
## where each .nc data file is a grid containing the 
## precipitation values for all grid points in a single 
## day, into individual time series; each new file 
## corresponds to the entire time series for a single 
## grid point. This script takes a long time to run 
## (around 1 day) and costs around 30GB of storage
############################################################

## Make sure user really wants to run this file
user_input <- readline(prompt = "Warning: This script will create 576,000 files. Continue? (Enter Y or N): ")
if(user_input != "Y"){
  stop("Program terminated")
}

###############################
##  Load necessary libraries ##
###############################

library(RNetCDF) ## For handling .nc files

###############################
##  Step 1: Import grid data ##
###############################

## First, set path to where this R file is
setwd("~/Research/DataConversion")

## Change directory to where .nc files are
## The nc files should be in a folder inside the folder
## containing this R file
main_dir   <- getwd()
setwd(paste(main_dir, "/netCDF_files", sep = ""))

## IMPORTANT -- There are too many files to load at once,
## so I'm looping through smaller chunks at a time - this
## is what my k index is doing. You can make k loop through
## a longer index to reduce the memory load (e.g. for(k in 1:6)),
## but you would need to adjust the index on i to reflect
## this change
filenames  <- list.files(pattern = "*.nc")
numDays    <- length(filenames)


for(k in 1:3){
  cat(k, "\n")
  setwd(paste(main_dir, "/netCDF_files", sep = ""))
  total.data <- list()
  #for(i in (1000 * (k - 1) + 1):min((k * 1000), numDays)){
  for(i in (2000 * (k - 1) + 1):min((k * 2000), numDays)){

    ## Open and read in file
    ncfile          <- filenames[i]
    nc.data         <- open.nc(ncfile, write = TRUE, path = "Data")
    cur.data        <- read.nc(nc.data, unpack = TRUE)
    
    ## Function to move all columns right; I'm converting 
    ## longitude from (0, 360) to (-180, 180)
    prec <- cur.data[[4]]
    prec_new <- matrix(0, nrow(prec), ncol(prec))
    prec_new[(((1:nrow(prec))+719) %% 1440) + 1, ] <- prec[1:nrow(prec),]
    
    ## Save data for the whole map
    total.data[[i - 2000 * (k - 1)]] <- prec_new
    close.nc(nc.data)
    
  }
  setwd(main_dir)
  
  ##################################################
  ##  Step 2: Convert into individual time series ##
  ##################################################
  
  main_dir   <- getwd()
  for(i in 1:1440){
    cat("Created folder", i, "\n")
    
    ## Create new folder for each longitude value
    new_folder <- paste("Long", i, sep = "")

    ## Create folder only for the first chunk (i.e. k = 1)
    if(k == 1){
      dir.create(new_folder, showWarnings = FALSE) 
    }
    
    ## Change directory to newly created folder
    setwd(paste(main_dir, "/", new_folder, sep = ""))
    
    ## For each longitude folder, save latitude values as separate files
    for(j in 1:400){
      indData   <- unlist(lapply(total.data, function(x) x[i, j]))
      fname     <- paste("Long_", i, "_Lat_", j, ".txt", sep = "")
      write.table(indData, file = fname, append = TRUE, 
                  row.names = FALSE, col.names = FALSE)
    }
    
    ## Go back to main directory
    setwd(main_dir)
    
  }

}
