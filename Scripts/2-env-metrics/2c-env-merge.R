################
# Mask out land and everything below 60 S lat for env var data and calc PCA
#
################


# Specify dir --------------------------------------------------

#path CY machine
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)


# read in data -------------------------------------------------

#environmental variability - add relative slope
env_var <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv')) %>%
  dplyr::mutate(rel_slope = slope / sd_resid)

#monthly spectral exponent
lf <- list.files(paste0(dir, 'data/L1/climate/era5/'), full.names = TRUE)
ef <- grep('spectral', lf, value = TRUE)
#read all files in and bind - change name to match env_var
se <- do.call(rbind, lapply(ef, read.csv)) %>%
  reshape2::melt(id.vars = c('cell_id', 'lon', 'lat'),
                 value.name = 'spectral_exp', variable.name = 'var')
se$var <- as.character(se$var)
se$var[which(se$var == 'spectral_beta_temp')] <- 'temp'
se$var[which(se$var == 'spectral_beta_precip')] <- 'precip'

#landmask
lm <- sf::st_read(paste0(dir, 'data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')


# apply land mask -------------------------------------------------------------

#just need one mask for both vars
frast <- dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, cell_id) %>%
  terra::rast(crs = "epsg:4326") %>%
  #mask out water
  terra::mask(terra::vect(lm))

#coords (lat/lon)
t_crds <- terra::crds(frast)
#values
t_vals <- data.frame(terra::values(frast)) %>%
  dplyr::filter(!is.na(cell_id))

#masked df (valid col)
#combine, filter above -60 lat, add valid col
t_cval <- cbind(t_crds, t_vals) %>%
  dplyr::rename(lat = y, lon = x) %>%
  dplyr::filter(lat > -60) %>%
  dplyr::mutate(valid = TRUE)


# merge env data ----------------------------------------------------------

#merge masked df with actual data
#set all ocean and land < 60 S lat to FALSE for field valid
env_mrg <- dplyr::left_join(env_var, se, 
                          by = c('cell_id', 'var', 'lon', 'lat')) %>%
  dplyr::left_join(t_cval, by = c('cell_id', 'lon', 'lat'))
na_idx <- which(is.na(env_mrg$valid))
env_mrg$valid[na_idx] <- FALSE

# saveRDS(env_mrg, '~/tt.rds')


# PCA ---------------------------------------------------------------------

#function to get correlation between metrics, conduct pca, and add pc to env data
pca_fun <- function(df, VAR)
{
  tenv <- dplyr::filter(df, valid == TRUE, var == VAR)
  tcor <- round(cor(cbind(tenv$sd_resid, tenv$sd_season, tenv$mean)), 3)
  
  tcov <- tenv[,c('sd_resid', 'sd_season', 'mean')]
  colnames(tcov) <- c('Inter-annual sd', 'Intra-annual sd', 'Mean')
  covs_pca <- prcomp(tcov, center = TRUE, scale. = TRUE)
  
  #extract PCs and merge with data
  tenv$PC1 <- covs_pca$x[,1]
  tenv$PC2 <- covs_pca$x[,2]
  tenv$PC3 <- covs_pca$x[,3]
  
  tlist <- list(cor = tcor, covs_pca = covs_pca, mrg_out = tenv)
  return(tlist)
}

#run function
pc_temp <- pca_fun(df = env_mrg, VAR = 'temp')
pc_precip <- pca_fun(df = env_mrg, VAR = 'precip')
pc_temp_GAM <- pca_fun(df = env_mrg_GAM, VAR = 'temp')
pc_precip_GAM <- pca_fun(df = env_mrg_GAM, VAR = 'precip')


#pca plots
# pdf(paste0(fig_dir,'pca_biplot-', run_date, '.pdf'),
#     height = 5.75, width = 5.75)
factoextra::fviz_pca_var(pc_temp$covs_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         #col.var = "contrib", # Color by contributions to the PC
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE,
                         title = 'PCA - Temperature')
# dev.off()

factoextra::fviz_pca_var(pc_temp$covs_pca,
                         axes = c(2,3),
                         repel = FALSE,
                         title = 'PCA - Temperature')

factoextra::fviz_pca_var(pc_precip$covs_pca,
                         axes = c(1,2),
                         repel = FALSE,
                         title = 'PCA - Precip')

factoextra::fviz_pca_var(pc_precip$covs_pca,
                         axes = c(2,3),
                         repel = FALSE,
                         title = 'PCA - Precip')

#merge both temp and precip into main dfs
env_main <- rbind(pc_temp$mrg_out, pc_precip$mrg_out)
env_main_GAM <- rbind(pc_temp_GAM$mrg_out, pc_precip_GAM$mrg_out)


# PCA for temp/precip jointly ---------------------------------------------

tenv <- dplyr::filter(env_mrg, valid == TRUE, var == 'precip')
tcor <- round(cor(cbind(tenv$sd_resid, tenv$sd_season, tenv$mean)), 3)

hist(sqrt(tenv$mean))
hist(tenv$sd_resid)
hist(tenv$sd_season)

tcov <- tenv[,c('sd_resid', 'sd_season', 'mean')]
colnames(tcov) <- c('Inter-annual sd', 'Intra-annual sd', 'Mean')
covs_pca <- prcomp(tcov, center = TRUE, scale. = TRUE)

#extract PCs and merge with data
tenv$PC1 <- covs_pca$x[,1]
tenv$PC2 <- covs_pca$x[,2]
tenv$PC3 <- covs_pca$x[,3]

tlist <- list(cor = tcor, covs_pca = covs_pca, mrg_out = tenv)


# write out merged data -----------------------------------------------------

write.csv(env_main, paste0(dir, 'Data/L2/climate/era5/Env-main.csv'))
write.csv(env_main_GAM, paste0(dir, 'Data/L2/climate/era5/Env-main-GAM.csv'))
