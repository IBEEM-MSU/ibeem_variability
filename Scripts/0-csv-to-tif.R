# TITLE:            CSV <-> Tif Conversion File
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATE:             April 2024 
# OVERVIEW:         File to convert csvs back into tifs for use in data analysis. 

# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

# load libraries ------------------------------------------------

library(terra)
library(tidyverse)

# load tif data sets to convert ------------------------------------------------
tif_to_csv <- function(full_dir, filename, save_output = TRUE){
  csv <- terra::rast(paste0(full_dir, filename, ".tif")) %>% 
    terra::as.data.frame(., xy=TRUE, na.rm=FALSE)
  if(save_output == TRUE){
    write.csv(csv, paste0(full_dir, filename, ".csv"), row.names=F)
  }
  return(csv)
}

csv_to_tif <- function(full_dir, filename, save_output = TRUE){
  ras <- read.csv(paste0(full_dir, filename, ".csv")) %>% 
         terra::rast(., type="xyz", crs="EPSG:4326")
  if(save_output == TRUE){
    terra::writeRaster(ras, paste0(full_dir, filename, ".tif")) 
  }
  return(ras)
}

# load tif data sets to convert to csvs ------------------------------------------------
# tif_to_csv(full_dir = paste0(dir, "data/L3/"), filename = "delta")
# tif_to_csv(full_dir = paste0(dir, "data/L3/"), filename = "env_var")

# load csv data sets to convert to tifs ------------------------------------------------
csv_to_tif(full_dir = paste0(dir, "data/L3/"), filename = "delta", save_output = TRUE)
csv_to_tif(full_dir = paste0(dir, "data/L3/"), filename = "env_var", save_output = TRUE)
