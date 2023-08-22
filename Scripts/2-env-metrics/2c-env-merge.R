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
  #NOTE this cell_id is not the same as env_var cell_id (so remove it)
  dplyr::select(-cell_id) %>%
  reshape2::melt(id.vars = c('lon', 'lat'),
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
  terra::mask(terra::vect(lm), touches = TRUE)

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
  #dplyr::filter(lat > -66.5, lat < 66.5) %>%
  dplyr::mutate(valid = TRUE)


# merge env data ----------------------------------------------------------

#merge masked df with actual data
#set all ocean and land < 60 S lat to FALSE for field valid
env_mrg <- dplyr::left_join(env_var, se, 
                          by = c('var', 'lon', 'lat')) %>%
  dplyr::left_join(t_cval, by = c('cell_id', 'lon', 'lat')) %>%
  dplyr::rename(sp_color_month = spectral_exp)
na_idx <- which(is.na(env_mrg$valid))
env_mrg$valid[na_idx] <- FALSE

# saveRDS(env_mrg, '~/tt.rds')

# plt <- dplyr::filter(env_mrg, var == 'temp', valid == TRUE) %>%
#   dplyr::select(lon, lat, sp_color_yearly) %>%
#   terra::rast(crs = "epsg:4326")
# plot(plt)


# explore -----------------------------------------------------------------

tt <- dplyr::filter(env_mrg, var == 'precip', valid == TRUE,
                    !is.na(sp_color_yearly))
# plot(tt$sd_resid, tt$sd_season, col = rgb(0,0,0,0.1))
# plot(tt$sd_resid, tt$sp_color_monthly, col = rgb(0,0,0,0.1))
# plot(tt$sd_resid, tt$sp_color_year, col = rgb(0,0,0,0.1))
# plot(tt$sd_season, tt$sp_color_monthly, col = rgb(0,0,0,0.1))
# plot(tt$sd_season, tt$sp_color_yearly, col = rgb(0,0,0,0.1))
# plot(tt$sp_color_yearly, tt$sp_color_monthly, col = rgb(0,0,0,0.1))

cor(cbind(tt$mean, tt$sd_resid, tt$sd_season, tt$sp_color_yearly, tt$sp_color_monthly))
tcov <- tt[,c('mean', 'sd_resid', 'sd_season', 'sp_color_yearly', 'sp_color_monthly')]
colnames(tcov) <- c('Mean', 'Inter-annual', 'Intra-annual', 'Spec (year)', 'Spec (month)')
covs_pca <- prcomp(tcov, center = TRUE, scale. = TRUE)

factoextra::fviz_pca_var(covs_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE,
                         title = 'PCA - Temperature')

tt$pc1 <- covs_pca$x[,1]
tt$pc2 <- covs_pca$x[,2]
tt$pc3 <- covs_pca$x[,3]

plt <- dplyr::select(tt, lon, lat, mean) %>%
  terra::rast(crs = "epsg:4326")
plot(plt)


#both temp and precip together
tt_temp <- dplyr::filter(env_mrg, var == 'temp')
names(tt_temp) <- c('cell_id', 'var', 'lon', 'lat', 
                    paste0('temp_', names(tt_temp)[-c(1:4)]))
tt_precip <- dplyr::filter(env_mrg, var == 'precip') %>%
  dplyr::select(-c(var, valid, cell_id, lon, lat))
names(tt_precip) <- paste0('precip_', names(tt_precip))
tt_mrg <- cbind(tt_temp, tt_precip) %>%
  dplyr::rename(valid = temp_valid)
str(tt_mrg)

tt_f <- dplyr::filter(tt_mrg, valid == TRUE, 
                      !is.na(precip_sp_color_yearly)) %>%
  dplyr::select(lon, lat, 
    # #TEMP
    temp_mean,
    # precip_mean,
    # #INTER-ANNUAL SD
    temp_sd_resid,
    # precip_sd_resid,
    # #INTRA-ANNUAL SD
    # temp_sd_season,
    # precip_sd_season,
    # #YEAR SPECTRAL COLOR
    temp_sp_color_yearly,
    # precip_sp_color_yearly,
    # #MONTH SPECTRAL COLOR
    # temp_sp_color_monthly,
    # precip_sp_color_monthly
    )

cor(tt_f)
covs_pca <- dplyr::select(tt_f, -lon, -lat) %>%
  prcomp(center = TRUE, scale. = TRUE)

factoextra::fviz_pca_var(covs_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE,
                         title = 'PCA - Temp and Precip')
tt_f$pc1 <- covs_pca$x[,1]
tt_f$pc2 <- covs_pca$x[,2]

plt <- dplyr::select(tt_f, lon, lat, pc1) %>%
  # dplyr::mutate(tt = scale(temp_sd_resid, scale = TRUE)) %>%
  # dplyr::mutate(tt = log(temp_sd_resid)) %>%
  dplyr::mutate(tt = pc1) %>%
  dplyr::select(lon, lat, tt) %>%
  terra::rast(crs = "epsg:4326")
rv <- range(terra::values(plt$tt), na.rm = TRUE)
pal <- leaflet::colorNumeric(palette = "RdBu", 
                             domain = rv, 
                             reverse = T)
plot(plt, range = rv, col = pal(seq(rv[1], rv[2], by = 0.1)))



# PCA ---------------------------------------------------------------------

#variables: mean, sd_year, sd_season, sp_color_year, sp_color_month

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

# t1 <- dplyr::filter(env_mrg, var == 'temp')
# str(t1)
# plot(t1$spectral_exp, t1$spectral_beta)
# range(t1$spectral_beta, na.rm = TRUE)
# 
# #make wide rather than long
# #also look at correlation betwee
# tenv <- dplyr::filter(env_mrg, valid == TRUE, var == 'precip')
# tcor <- round(cor(cbind(tenv$sd_resid, tenv$sd_season, tenv$mean)), 3)
# 
# hist(sqrt(tenv$mean))
# hist(tenv$sd_resid)
# hist(tenv$sd_season)
# 
# tcov <- tenv[,c('sd_resid', 'sd_season', 'mean')]
# colnames(tcov) <- c('Inter-annual sd', 'Intra-annual sd', 'Mean')
# covs_pca <- prcomp(tcov, center = TRUE, scale. = TRUE)
# 
# #extract PCs and merge with data
# tenv$PC1 <- covs_pca$x[,1]
# tenv$PC2 <- covs_pca$x[,2]
# tenv$PC3 <- covs_pca$x[,3]
# 
# tlist <- list(cor = tcor, covs_pca = covs_pca, mrg_out = tenv)


# write out merged data -----------------------------------------------------

write.csv(env_main, paste0(dir, 'Data/L2/climate/era5/Env-main.csv'))
