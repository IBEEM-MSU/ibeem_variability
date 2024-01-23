####################
# Figure 5 - delta metrics
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
#XXXX


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


