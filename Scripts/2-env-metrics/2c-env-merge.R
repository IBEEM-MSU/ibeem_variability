################
# Merge seasonality and env
#
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)


# read in data -------------------------------------------------

#environmental variability
env_var <- read.csv(paste0(dir, 'Data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))
env_var_GAM <- read.csv(paste0(dir, 'Data/L2/climate/era5/Env-var-GAM-1_2_3_4_5_6_7_8_9_10_11_12.csv'))
#environmental seasonality
env_season <- read.csv(paste0(dir, 'Data/L1/climate/era5/Env-seasonality-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

#landmask
lm <- sf::st_read(paste0(dir, 'Data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')


# merge env data ----------------------------------------------------------

#merge
temp_mrg <- dplyr::filter(env_var, var == 'temp') %>%
  dplyr::left_join(env_season, by = c('lat', 'lon')) %>%
  dplyr::select(-mn_sd_precip) %>%
  dplyr::rename(sd_season = mn_sd_temp)

precip_mrg <- dplyr::filter(env_var, var == 'precip') %>%
  dplyr::left_join(env_season, by = c('lat', 'lon')) %>%
  dplyr::select(-mn_sd_temp) %>%
  dplyr::rename(sd_season = mn_sd_precip)

env_mrg <- rbind(temp_mrg, precip_mrg)


# rasterize and mask env var data ----------------------------------------------

#multiband raster (each band different env metric)
#var = temp or precip
#function
mbr_fun <- function(input, VAR)
{
  #rasterize
  trast <- dplyr::filter(input, var == VAR) %>%
    dplyr::select(-cell_id, -var) %>%
    terra::rast(crs = "epsg:4326") %>%
    #mask out water
    terra::mask(terra::vect(lm))
  
  #get relative slope
  sl_resid <- trast[['slope']] / trast[['sd_resid']]
  names(sl_resid) <- 'rel_slope'
  
  #add band
  trast2 <- c(trast, sl_resid)
  
  return(trast2)
}


#run function - putting some cells outside of lat/lon bounds for some reason
ev_temp <- mbr_fun(input = env_mrg, VAR = 'temp')
ev_precip <- mbr_fun(input = env_mrg, VAR = 'precip')
# ev_temp_GAM <- mbr_fun(input = env_var_GAM, VAR = 'temp')
# ev_precip_GAM <- mbr_fun(input = env_var_GAM, VAR = 'precip')

#bands
names(ev_temp)

#check
# plot(ev_temp[[1]])
# plot(ev_temp[[4]])
# plot(ev_temp[[14]])


# extract values from masked raster ---------------------------------------

#coords (lat/lon)
t_crds <- terra::crds(ev_temp)
#values
t_vals <- data.frame(terra::values(ev_temp)) %>%
  dplyr::filter(!is.na(mean))
#combine
t_mrg <- cbind(t_crds, t_vals) %>%
  dplyr::rename(lat = y, lon = x) %>%
  dplyr::filter(lat > -60)


# PCA ---------------------------------------------------------------------

cor(cbind(t_mrg$sd_resid, t_mrg$mean, t_mrg$sd_season))

t_cov <- t_mrg[,c('sd_resid', 'mean', 'sd_season')]
colnames(t_cov) <- c('Inter-annual sd', 'Mean', 'Intra-annual sd')
covs_pca <- prcomp(t_cov, center = TRUE, scale. = TRUE)

# pdf(paste0(fig_dir,'pca_biplot-', run_date, '.pdf'),
#     height = 5.75, width = 5.75)
factoextra::fviz_pca_var(covs_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         #col.var = "contrib", # Color by contributions to the PC
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE)
# dev.off()

factoextra::fviz_pca_var(covs_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         #col.var = "contrib", # Color by contributions to the PC
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE)

#extract PCs
t_mrg$PC1 <- covs_pca$x[,1]
t_mrg$PC2 <- covs_pca$x[,2]
t_mrg$PC3 <- covs_pca$x[,3]


# rerasterize -------------------------------------------------------------

tt <- terra::rast(t_mrg, crs = "epsg:4326")
names(tt)

#mean
plot(tt[[1]], main = 'mean')
plot(tt[[4]], main = 'Inter-annual sd')
plot(tt[[14]], main = 'Intra-annual sd')

plot(tt[[16]], main = 'PC1 (+ = ^ mean, v sd)')
plot(tt[[17]], main = 'PC2 (+ = ^ intra, v inter)')
plot(tt[[18]], main = 'PC3 (+ = ^ mean, ^ intra)')


# stats -------------------------------------------------------------------

#median over globe
terra::global(ev_temp[['slope']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['sd_resid']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['kurt']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['skew']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['spectral_beta']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['rho_l1']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['rel_slope']], fun = function(x) median(x, na.rm = TRUE))


# plots -----------------------------------------------------------------

# NA_out <- dplyr::filter(env_var, 
#                            lon > -170, lon < -50,
#                            lat > 15, lat < 75)

#function to plot lm and gam
gc_fun <- function(rast, rast_gam, var)
{
  par(mfrow = c(2,1))
  
  #get range of values  
  rng <- range(terra::global(c(rast[[var]], rast_gam[[var]]), 
                             fun = function(x) range(x, na.rm = TRUE)))
  plot(rast[[var]], 
       main = var,
       range = rng)
  plot(rast_gam[[var]], 
       main = paste0(var, ' - GAM'),
       range = rng)
}

#run fun
#slope doesn't vary since calcated with lm for both
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'slope')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'sd_resid')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'kurt')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'skew')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'spectral_beta')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'rho_l1')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'rel_slope')


#which areas are highly predictable (P; high temporal autocorrelation) on short time scales (S) and have low intrinsic variability (IV)?
#high P, short S, low IV = short LH
#low P, long S, high IV = long LH

#high P/short S = faster LH
#high IV = slower LH

#is P at lag = gen time the same across all species?

