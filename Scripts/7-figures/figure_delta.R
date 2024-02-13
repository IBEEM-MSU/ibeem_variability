####################
# delta metrics (histograms and maps)
####################


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# dir <- '/mnt/research/ibeem/variability/'
# sc_dir <- '/mnt/home/ccy/variability/'
gl_run_date <- '2023-10-17'


# read in data ------------------------------------------------------------

#df from results
bird_df <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
                          '/bird-gl-phylo-vint-data-', gl_run_date, '.rds'))$pro_data

#stacked raster
del_ras <- terra::rast(paste0(dir, 'data/L3/raster-gl-dT-dP-nsp.tif')) 


# delta histograms ------------------------------------------------------

dt_hist <- hist(bird_df$temp_delta, breaks = 50, plot = F)
dt_col <- ifelse(dt_hist$breaks >= 0.1 & dt_hist$breaks < 0.3, rgb(0.2, 0.8, 0.5, 0.5), 
                 ifelse(dt_hist$breaks >= 0.3, 'purple', rgb(0.2,0.2,0.2,0.2)))
# range(bird_df$temp_delta)
# XLIM <- c(-0.2, 1)
plot(dt_hist, col = dt_col, border = FALSE, main = 'delta_T', 
     xlab = 'Rate of change degrees C (sd / generation)')#, 
#      xlim = XLIM,
#      xaxt = 'n')
# axis(1, at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))

dp_hist <- hist(bird_df$precip_delta, breaks = 50, plot = F)
dp_col <- ifelse((dp_hist$breaks >= 0.1 & dp_hist$breaks < 0.3) | dp_hist$breaks <= -0.1, rgb(0.2, 0.8, 0.5, 0.5), 
                 ifelse(dp_hist$breaks >= 0.3, 'purple', rgb(0.2,0.2,0.2,0.2)))
plot(dp_hist, col = dp_col, border = FALSE, main = 'delta_P', 
     xlab = 'Rate of change mm precipitation (sd / generation)')#,
# xlim = XLIM,
# xaxt = 'n')
# range(bird_df$precip_delta)
# axis(1, at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))


# delta maps --------------------------------------------------------------

# Mask areas with fewer than NSP species
NSP <- 10
msk <- ifel(del_ras[['n_sp']] < NSP, NA, 1)
del_ras2 <- terra::mask(del_ras, msk, 
                         inverse = FALSE)

# Land outline
# download land outline (50 m) from natural earth
land_50m <- ne_download(scale = 50,
                        category = "physical",
                        type = "land",
                        returnclass = c("sf"))

# Map bounding box

bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), 
          crs = st_crs(4326)),
  n = 100))

### Transform data to Robinson projection ----

# Transform raster to Robinson ("ESRI:54030")
del_ras_proj <- terra::project(del_ras2, "ESRI:54030")

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs = "ESRI:54030", type = "robin")

### Median gl map ----

# Make data frame with projected env data 
del_df <- terra::as.data.frame(del_ras_proj, xy = T, cells = F, na.rm = T)
head(del_df)

med_dT_q <- quantile(del_df$median_dT, seq(0, 1, by = 0.05))
med_dP_q <- quantile(del_df$median_dP, seq(0, 1, by = 0.05))
n_sp <- quantile(del_df$n_sp, seq(0, 1, by = 0.05))

# Map with ggplot 
base <- ggplot() +
  # Map background
  geom_sf(data = land_50m_robinson,
          color = "gray80",
          fill = "gray95",
          linetype = "solid",
          size = 0.2) +
  geom_sf(data = bb_robinson,
          color = "white",
          fill = NA,
          linetype = "solid") +
  theme_minimal()
  

# delta T
base + geom_tile(data = del_df, 
            aes(x = x,
                y = y, 
                fill = median_dT)) +
  scale_fill_gradient(low = "#fef0d9",
                      high = "#b30000",
                      na.value = "#FAFAFA",
                      #limits = c(0,5)) +
                      limits = range(med_dT_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle('Delta T')

# delta P
base + geom_tile(data = del_df, 
                 aes(x = x,
                     y = y, 
                     fill = median_dP)) +
  scale_fill_gradient(low = "#fef0d9",
                      high = "#b30000",
                      na.value = "#FAFAFA",
                      #limits = c(0,5)) +
                      limits = range(med_dP_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle('Delta P')


# correlation delta T and delta P -----------------------------------------

cor(bird_df$precip_delta, bird_df$temp_delta)
