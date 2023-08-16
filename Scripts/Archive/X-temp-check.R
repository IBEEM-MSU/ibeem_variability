
# Specify dir --------------------------------------------------

dir <- '~/Downloads/era5/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)


# read in data -------------------------------------------------

#read in env data
env_var <- read.csv(paste0(dir, 'Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

#read in climate data
cl <- read.csv('~/Desktop/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv')


# process -----------------------------------------------------------------

#hist spectral exponent values
tt <- dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, spectral_beta)
hist(tt$spectral_beta)

#filter for cellls with negative spectral exponents
#add cell id
nse <- dplyr::filter(env_var, var == 'temp', spectral_beta < 0) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup()

udf <- unique(nse[,c('lat', 'lon', 'cell_id')])

#filter climate data for those cells
fcl <- dplyr::left_join(cl, udf) %>%
  dplyr::filter(!is.na(cell_id))

#unique cells
uci <- unique(fcl$cell_id)

#look at those time series
#i <- 1
t_fcl2 <- dplyr::filter(fcl2, cell_id == uci[i]) %>%
  arrange(year)
plot(t_fcl2$year, t_fcl2$temp, type = 'l')

te <- t_fcl2
#linear model fit
fit_temp <- summary(lm(temp ~ year, data = te))

#residuals from model
temp_resid <- residuals(fit_temp)
plot(temp_resid, type = 'l')
plot(te$year, temp_resid, type = 'l')

#spectral analysis using Lomb-Scargle Periodogram
#between freq 2/(n*dt) and 1/(2*dt), where dt = 1 and n = 72; 0.0278 to 0.5
#following Marshall and Burgess 2015 Eco Letters
#see also Vasseur and Yodzis 2004 Ecology
#period = 1/freq; ~2 - 36 years
temp_spec <- lomb::lsp(temp_resid, from = 0.0278, to = 0.5, type = 'frequency',
                       normalize =  'standard', plot = FALSE)
plot(temp_spec)
#spectral exponent (1/f^beta)
(temp_spec_fit <- summary(lm(log10(temp_spec$power) ~ log10(temp_spec$scanned)))$coefficients[,1])

#temporal autocorrelation
temp_acf <- acf(temp_resid, lag.max = 5, plot = FALSE)
precip_acf <- acf(precip_resid, lag.max = 5, plot = FALSE)



