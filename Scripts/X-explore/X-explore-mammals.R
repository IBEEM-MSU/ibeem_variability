################
# Initial exploration mammals
#
#
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(viridis)
library(lme4)


# read in data -------------------------------------------------

mam_df <- read.csv(paste0(dir, 'Data/L3/main-mammal-data.csv')) %>%
  dplyr::mutate(Family = PH_Family,
                Order = PH_Order,
                # LH_Mass = LH_AdultBodyMass_g, #Pacifici mass
                Mass = PH_Mass.g, #Phylacine mass
                GenLength = LH_GenerationLength_d,
                fac_Family = factor(Family),
                fac_Order = factor(Order),
                lMass = log(Mass),
                lGL = log(GenLength))

#one species with Inf precip_cv_season (due to precip_mean = 0 and using median)
dplyr::filter(mam_df, !is.finite(precip_cv_season))
mam_df$precip_cv_season[which(!is.finite(mam_df$precip_cv_season))] <- 0


# raw temp and precip -----------------------------------------------------

#candidate models
f1 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year +
                   temp_sp_color_month +
                   precip_cv_season +
                   precip_cv_year +
                   precip_sp_color_month +
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sp_color_month | fac_Family) +
                   (-1 + precip_cv_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family) +
                   (-1 + precip_sp_color_month | fac_Family),
                 data = mam_df)
f2 <- lme4::lmer(lGL ~ lMass + 
                   precip_cv_season + 
                   precip_cv_year +
                   (1 + lMass | fac_Family) +
                   (-1 + precip_cv_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = mam_df)
f3 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year +
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = mam_df)
f4 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = mam_df)
f5 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   precip_cv_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = mam_df)
f6 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = mam_df)
f7 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = mam_df)
round(summary(f1)$coefficients, 2)

fx <- lm(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year +
                   # temp_sp_color_month +
                   precip_cv_season +
                   precip_cv_year, #+
                   # precip_sp_color_month,
                 data = mam_df)
summary(fx)
car::vif(fx)

AIC(f1) #2nd best
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6) #best ()
AIC(f7)

summary(f1)
summary(f4)
summary(f5)
summary(f6)

#VIF - by any metric, looks good (VIF < 5 is OK; < 3 stricter threshold?)
car::vif(f1)
