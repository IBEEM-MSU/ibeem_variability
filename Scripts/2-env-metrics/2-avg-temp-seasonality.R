#---------------------------------------------------------
# Calculate average temperature and seasonality per cell

# Load packages -----------------------------

library(tidyverse)

# Read in data ------------------------------
clim_dir_in <- '/mnt/research/ibeem/data/L1/climate/era5/'
clim_dir_out <- '/mnt/research/ibeem/L2/climate/era5/'
env_out <- read.csv(paste0(clim_dir_in, 'ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

# Calculate average temp --------------------
#add unique cell id
#convert to data table (much faster)
env_out2 <- dplyr::group_by(env_out, lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cell_id, year) %>%
  data.table::as.data.table()

# Calculate average temp per cell
avg_temp <- env_out2 %>%
  dplyr::group_by(cell_id) %>%
  dplyr::reframe(lat = as.numeric(lat),
          lon = as.numeric(lon),
          avg_temp = mean(temp))

# Save file ----------------------------------
write.csv(avg_temp, file= paste0(clim_dir_out, 'avg_temp.csv'), row.names = FALSE)
