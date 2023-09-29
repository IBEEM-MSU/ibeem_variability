################
# Initial exploration bird
#
#
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(viridis)
library(lme4)


# read in data -------------------------------------------------

#bird df
or_excl <- c('Sphenisciformes', #penguins 
             'Procellariiformes', #tubenoses
             'Pelecaniformes', #pelicans
             'Suliformes', #gannets/boobies
             'Phaethontiformes', #tropicbirds
             'Charadriiformes')#, #skuas, gulls, terns, skimmers, auks
#'Anseriformes', #waterfowl
#'Ciconiiformes', #storks
#'Gaviiformes', #aquatic birds (loons and divers)
#'Gruiformes', #cranes, rails - Family Psophiidae are not waterbirds, but there are few members (~6 species)
#'Phoenicopteriformes', #flamingos and relatives
#'Podicipediformes') #grebes

'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1) %>%
  dplyr::mutate(fac_Family = factor(Family),
                fac_Order = factor(Order),
                lMass = log(Mass),
                lGL = log(GenLength))

#bird raster
bird_ras <- terra::rast(paste0(dir, 'data/L3/raster-gl-dh-nsp.tif'))
bird_ras2 <- bird_ras[[c('median_gl', 'sd_gl', 
                   'median_dh', 'sd_dh')]]
#mask areas with fewer than 5 species
msk <- ifel(bird_ras[['n_sp']] < 5, NA, 1)
bird_ras3 <- terra::mask(bird_ras2, msk, 
                         inverse = FALSE)


# map gen length (species median lat/lon) ---------------------------------------------------

#USE RASTER INSTEAD
# #map
# world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# 
# #base plot
# base_plt <- ggplot() + 
#   #plot map
#   geom_sf(data = world,
#           #fill = 'grey89') +
#           fill = 'aliceblue') +
#   # theme_bw()
#   #graticles
#   theme(panel.grid.major = #element_line(color = black,
#           #            linetype = 'dashed', 
#           #           linewidth = 0.1),
#           element_blank(),
#         panel.ontop = TRUE,
#         panel.background = element_rect(fill = 'NA')) +
#   scale_x_continuous(breaks = seq(-180, 180, by = 15))+
#   scale_y_continuous(breaks = seq(-80, 80, by = 20))
# 
# #pts
# gen_pts_wo_seabirds_coastal <- base_plt +
#   #plot pts
#   geom_point(data = bird_df,
#              aes(cen_lon, cen_lat,
#                  col = lGL),
#              alpha = 0.5,
#              size = 1.1) +
#   scale_color_viridis()
# 
# gen_pts_all <- base_plt +
#   #plot pts
#   geom_point(data = bird_df,
#              aes(cen_lon, cen_lat,
#                  col = lGL),
#              alpha = 0.5,
#              size = 1.1) +
#   scale_color_viridis()
# 
# gen_pts_wo_seabirds_coastal
# gen_pts_all


# map gen length raster ---------------------------------------------------

plot(bird_ras3[[c('median_gl', 'median_dh')]])
plot(bird_ras3[[c('sd_gl', 'sd_dh')]])
plot(bird_ras[['n_sp']])


# how long before temp trend exceeds 2 sd interannual ---------------------

plot(bird_df$temp_sd_year, bird_df$temp_slope, 
     col = rgb(0,0,0,0.1),
     pch = 19)

#ex = number of years before temp exceeds 2 sd
#AND
#delta_t = how much temp will change in 1 generation (in sds)
#AND
#delta_haldane = how much change (in sd) per generation
#AND
#n_gen = how many gens before temp will exceed 2 sd
ex_df <- dplyr::mutate(bird_df, 
                       ex = 2 * temp_sd_year / temp_slope,
                       ex_season = 2 * temp_sd_season / temp_slope,
                       delta_t = temp_slope / temp_sd_year * GenLength,
                       delta_t_season = temp_slope / temp_sd_year * GenLength,
                       #(degrees / year) * (sd / degrees) * (year / gen) = sd / gen
                       delta_haldane = (temp_slope / temp_sd_year) * GenLength,
                       delta_haldane_season = (temp_slope / temp_sd_season) * GenLength,
                       n_gen = ex / GenLength)
hist(ex_df$ex, xlim = c(0, 125), breaks = 10000)
hist(ex_df$ex_season, xlim = c(0, 125), breaks = 100000)
plot(ex_df$ex, ex_df$GenLength, xlim = c(0, 125),
     col = rgb(0,0,0,0.1),
     pch = 19)
hist(ex_df$delta_t, breaks = 100)
hist(ex_df$n_gen, breaks = 10000, xlim = c(0, 125))
hist(log(ex_df$n_gen))
#sd's per generation
hist(ex_df$delta_haldane, 
     xlim = c(0, 0.5), breaks = 100,
     main = 'delta temp (sd / generation)')
abline(v = 0.1, lty = 2, col = 'red', lwd = 3)
abline(v = 0.3, lty = 2, col = 'red', lwd = 3)
hist(ex_df$delta_haldane_season, 
     xlim = c(0, 0.5), breaks = 100,
     main = 'delta temp season (sd / generation)')
abline(v = 0.1, lty = 2, col = 'red', lwd = 3)
abline(v = 0.3, lty = 2, col = 'red', lwd = 3)


#Bürger and Lynch 1995: rate of phenotypic evolution is considerably less than 10% per gen
#Gingerich 2009: rate of phenotypic evolution is 0.1 to 0.3 SD per generation


# raw temp -----------------------------------------------------------

#candidate models
f1 <- lme4::lmer(lGL ~ lMass+ 
                   temp_sp_color_month + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sp_color_month | fac_Family) + 
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
f2 <- lme4::lmer(lGL ~ lMass + 
                   temp_sp_color_year + 
                   temp_sd_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sp_color_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family), 
                 data = bird_df)
f3 <- lme4::lmer(lGL ~ lMass + 
                   temp_sp_color_year + 
                   temp_sd_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sp_color_year | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
f4 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family), 
                 data = bird_df)
f5 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
f6 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
f7 <- lme4::lmer(lGL ~ lMass + 
                   (1 + lMass | fac_Family), 
                 data = bird_df)

fx <- lm(lGL ~ lMass + temp_sp_color_month + temp_sd_season + temp_sd_year, 
         data = bird_df)

#model selection
AIC(f1) #best (color)
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6) #close to best (no color)
AIC(f7)
AIC(fx)
MuMIn::r.squaredGLMM(f6)

summary(f1)
summary(f4)
summary(f5)
summary(f6)
summary(fx)

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
# % increase in Gen Length for 1 unit change in sd season
(exp(summary(f6)$coefficient[3,1]) - 1) * 100
# % increase in Gen Length for 1 unit change in sd year
(exp(summary(f6)$coefficient[4,1]) - 1) * 100
# %increase in Gen Length over observed range of sd season
(exp(summary(f6)$coefficient[3,1] * 
       diff(range(bird_df$temp_sd_season))) - 1) * 100
# %increase in Gen Length over observed range of sd year
(exp(summary(f6)$coefficient[4,1] * 
       diff(range(bird_df$temp_sd_year))) - 1) * 100


#VIF - by any metric, looks good (VIF < 5 is OK; < 3 stricter threshold?)
car::vif(f1)
vif_1 <- lme4::lmer(lMass ~
                      temp_sd_season + 
                      temp_sd_year + 
                      (1 + temp_sd_season | fac_Family) +
                      (-1 + temp_sd_year | fac_Family), 
                    data = bird_df)
vif_2 <- lme4::lmer(temp_sd_season ~
                      lMass + 
                      temp_sd_year + 
                      (1 + lMass | fac_Family) +
                      (-1 + temp_sd_year | fac_Family), 
                    data = bird_df)
vif_3 <- lme4::lmer(temp_sd_year ~
                      lMass +
                      temp_sd_season + 
                      (1 + lMass | fac_Family) +
                      (-1 + temp_sd_season | fac_Family), 
                    data = bird_df)

vif_1 <- summary(lm(lMass ~ temp_sd_year + temp_sd_season, data = bird_df))
vif_2 <- summary(lm(temp_sd_year ~ lMass + temp_sd_season, data = bird_df))
vif_3 <- summary(lm(temp_sd_season ~ lMass + temp_sd_year, data = bird_df))
vif_4 <- summary(lm(temp_sd_year ~ temp_sd_season, data = bird_df))
1 / (1 - vif_1$adj.r.squared)
1 / (1 - vif_2$adj.r.squared)
1 / (1 - vif_3$adj.r.squared)
1 / (1 - vif_4$adj.r.squared)
1 / (1 - MuMIn::r.squaredGLMM(vif_1)[,1])
1 / (1 - MuMIn::r.squaredGLMM(vif_2)[,1])
1 / (1 - MuMIn::r.squaredGLMM(vif_3)[,1])

ufam <- unique(bird_df$Family)
tdf <- data.frame(Family = ufam,
           VIF1 = NA,
           VIF2 = NA,
           VIF3 = NA,
           VIF4 = NA,
           VIF5 = NA,
           VIF6 = NA,
           VIF7 = NA)
for (i in 1:length(ufam))
{
  #i <- 1
  tt <- dplyr::filter(bird_df, Family == ufam[i])
  
  # tt <- bird_df
  if (NROW(tt) >= 10)
  {
    tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season + temp_sp_color_month +
       precip_cv_year + precip_cv_season + precip_sp_color_month, data = tt)
    stf1 <- summary(tf1)
    tdf$VIF1[i] <- 1 / (1 - stf1$r.squared)
    VIF1a <- 1 / (1 - stf1$adj.r.squared)
    
    tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season + temp_sp_color_month +
                precip_cv_year + precip_cv_season + precip_sp_color_month, data = tt)
    stf2 <- summary(tf2)
    tdf$VIF2[i] <- 1 / (1 - stf2$r.squared)
    VIF2a <- 1 / (1 - stf2$adj.r.squared)
    
    tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year + #temp_sp_color_month +
                precip_cv_year + precip_cv_season + precip_sp_color_month, data = tt)
    stf3 <- summary(tf3)
    tdf$VIF3[i] <- 1 / (1 - stf3$r.squared)
    VIF3a <- 1 / (1 - stf3$adj.r.squared)
    
    tf4 <- lm(temp_sp_color_month ~ lMass + temp_sd_year + temp_sd_season +
                precip_cv_year + precip_cv_season + precip_sp_color_month, data = tt)
    stf4 <- summary(tf4)
    tdf$VIF4[i] <- 1 / (1 - stf4$r.squared)
    VIF4a <- 1 / (1 - stf4$adj.r.squared)
    
    tf5 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
                temp_sp_color_month + 
                precip_cv_season + precip_sp_color_month, 
              data = tt)
    stf5 <- summary(tf5)
    tdf$VIF5[i] <- 1 / (1 - stf5$r.squared)
    VIF5a <- 1 / (1 - stf5$adj.r.squared)
    
    tf6 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
                temp_sp_color_month + 
                precip_cv_year + precip_sp_color_month, 
              data = tt)
    stf6 <- summary(tf6)
    tdf$VIF6[i] <- 1 / (1 - stf6$r.squared)
    VIF6a <- 1 / (1 - stf6$adj.r.squared)
    
    tf7 <- lm(precip_sp_color_month ~ lMass + temp_sd_year + temp_sd_season + 
                temp_sp_color_month + 
                precip_cv_year + precip_cv_season, 
              data = tt)
    stf7 <- summary(tf7)
    tdf$VIF7[i] <- 1 / (1 - stf7$r.squared)
    VIF7a <- 1 / (1 - stf7$adj.r.squared)
  }
}

median(tdf$VIF1, na.rm = TRUE)
median(tdf$VIF2, na.rm = TRUE)
median(tdf$VIF3, na.rm = TRUE)
median(tdf$VIF4, na.rm = TRUE)
median(tdf$VIF5, na.rm = TRUE)
median(tdf$VIF6, na.rm = TRUE)
median(tdf$VIF7, na.rm = TRUE)


#partial residual plots function
library(remef)
pr_fun <- function(model, data, var_interest, var_names)
{
  #get partial resids
  yp <- remef::remef(model, 
                     fix = var_names[-grep(var_interest, var_names)],
                     ran = "all")
  #get coefs
  coef_fit <- coef(summary(model))
  #plot resids
  plot(bird_df[,grep(var_interest, colnames(bird_df))], 
       yp, 
       col = rgb(0,0,0,0.1), 
       pch = 19,
       main = var_interest,
       xlab = var_interest,
       ylab = 'Partial resid; log(GenLength)')
  #plot model fit
  abline(a = coef_fit[1,1], 
         b = coef_fit[grep(var_interest, row.names(coef_fit)), 1], 
         col = 'red',
         lwd = 2)
}


#Mass
pr_fun(model = f1, 
       var_interest = 'lMass',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year'))
#temp season
pr_fun(model = f1, 
       var_interest = 'temp_sd_season',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year'))
#temp inter annual
pr_fun(model = f1, 
       var_interest = 'temp_sd_year',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year'))


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
                         data = bird_df)
f2 <- lme4::lmer(lGL ~ lMass + 
                   precip_cv_season + 
                   precip_cv_year +
                   (1 + lMass | fac_Family) +
                   (-1 + precip_cv_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df)
f3 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year +
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
f4 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = bird_df)
f5 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   precip_cv_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df)
f6 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_year + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df)
f7 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = bird_df)
flm <- lm(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year +
                   temp_sp_color_month +
                   precip_cv_season + 
                   precip_cv_year +
                   precip_sp_color_month,
                 data = bird_df)
round(summary(f1)$coefficients, 2)

AIC(f1) #best (sp color temp is within 3)
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6)
AIC(f7)
AIC(flm)

summary(f1)
summary(f4)
summary(f5)
summary(f6)
summary(flm)


hist(coef(f1)$fac_Family$temp_sd_season)
hist(coef(f1)$fac_Family$temp_sd_year)
hist(coef(f1)$fac_Family$precip_cv_season)
hist(coef(f1)$fac_Family$precip_cv_year)

# tplt_end <- data.frame(param = c('lMass', 'temp_sd_season', 'temp_sd_year', 
#                                  'precip_cv_season', 'precip_cv_year'),
#                        # med = summary(f1)$coef[-1,1],
#                        # lci = summary(f1)$coef[-1,1] - 1.96 * summary(f1)$coef[-1,2],
#                        # uci = summary(f1)$coef[-1,1] + 1.96 * summary(f1)$coef[-1,2],
#                        med = c((exp(summary(f1)$coef[2,1] * 
#                                    diff(range(bird_df$lMass))) - 1) * 100,
#                                (exp(summary(f1)$coef[3,1] * 
#                                       diff(range(bird_df$temp_sd_season))) - 1) * 100,
#                                (exp(summary(f1)$coef[4,1] * 
#                                       diff(range(bird_df$temp_sd_year))) - 1) * 100,
#                                (exp(summary(f1)$coef[5,1] * 
#                                       diff(range(bird_df$precip_cv_season))) - 1) * 100,
#                                (exp(summary(f1)$coef[6,1] * 
#                                       diff(range(bird_df$precip_cv_year))) - 1) * 100),
#                        lci = c((exp((summary(f1)$coef[2,1] - 1.96 * summary(f1)$coef[2,2]) * 
#                                 diff(range(bird_df$lMass))) - 1) * 100,
#                          (exp((summary(f1)$coef[3,1] - 1.96 * summary(f1)$coef[3,2]) *
#                                 diff(range(bird_df$temp_sd_season))) - 1) * 100,
#                          (exp((summary(f1)$coef[4,1] - 1.96 * summary(f1)$coef[4,2]) *
#                                 diff(range(bird_df$temp_sd_year))) - 1) * 100,
#                          (exp((summary(f1)$coef[5,1] - 1.96 * summary(f1)$coef[5,2]) *
#                                 diff(range(bird_df$precip_cv_season))) - 1) * 100,
#                          (exp((summary(f1)$coef[6,1] - 1.96 * summary(f1)$coef[6,2]) *
#                                 diff(range(bird_df$precip_cv_year))) - 1) * 100),
#                        uci = c((exp((summary(f1)$coef[2,1] + 1.96 * summary(f1)$coef[2,2]) * 
#                                       diff(range(bird_df$lMass))) - 1) * 100,
#                                (exp((summary(f1)$coef[3,1] + 1.96 * summary(f1)$coef[3,2]) *
#                                       diff(range(bird_df$temp_sd_season))) - 1) * 100,
#                                (exp((summary(f1)$coef[4,1] + 1.96 * summary(f1)$coef[4,2]) *
#                                       diff(range(bird_df$temp_sd_year))) - 1) * 100,
#                                (exp((summary(f1)$coef[5,1] + 1.96 * summary(f1)$coef[5,2]) *
#                                       diff(range(bird_df$precip_cv_season))) - 1) * 100,
#                                (exp((summary(f1)$coef[6,1] + 1.96 * summary(f1)$coef[6,2]) *
#                                       diff(range(bird_df$precip_cv_year))) - 1) * 100))
# 
# 
# tplt_ss_end <- data.frame(lmass = tplt_end[1,2] + 
#                             (exp(coef(f1)$fac_Family$lMass) - 1) * 100,
#                           temp_sd_season = tplt_end[2,2] + 
#                             (exp(coef(f1)$fac_Family$temp_sd_season) - 1) * 100,
#                           temp_sd_year = tplt_end[3,2] + 
#                             (exp(coef(f1)$fac_Family$temp_sd_year) - 1) * 100,
#                           precip_cv_season = tplt_end[4,2] + 
#                             (exp(coef(f1)$fac_Family$precip_cv_season) - 1) * 100,
#                           precip_cv_year = tplt_end[5,2] + 
#                             (exp(coef(f1)$fac_Family$precip_cv_year) - 1) * 100)
# 
# pnt_plt <- ggplot(tplt_ss_end, aes(-1.5, temp_sd_season)) +
#   geom_jitter(shape = 1, alpha = 0.2, size = 3) +
#   geom_hline(yintercept = 0,
#              linetype = 'dashed',
#              linewidth = 1,
#              alpha = 0.4) +
#   geom_segment(aes(y = tplt_end$lci[2], yend = tplt_end$uci[2],
#                    x = -1.5, xend = -1.5), 
#                lineend = "round",
#                linewidth = 3, alpha = 0.025) +
#   geom_point(aes(-1.5, tplt_end$med[2]), 
#              size = 8, alpha = 0.025) +
#   #temp sd year
#   geom_jitter(data = tplt_ss_end, aes(-0.5, temp_sd_year),
#               shape = 1, alpha = 0.2, size = 3) +
#   geom_segment(aes(y = tplt_end$lci[3], yend = tplt_end$uci[3],
#                    x = -0.5, xend = -0.5), 
#                lineend = "round",
#                linewidth = 3, alpha = 0.025) +
#   geom_point(aes(-0.5, tplt_end$med[3]), 
#              size = 8, alpha = 0.025) +
#   #precip_cv_season
#   geom_jitter(data = tplt_ss_end, aes(0.5, precip_cv_season),
#               shape = 1, alpha = 0.2, size = 3) +
#   geom_segment(aes(y = tplt_end$lci[4], yend = tplt_end$uci[4],
#                    x = 0.5, xend = 0.5), 
#                lineend = "round",
#                linewidth = 3, alpha = 0.025) +
#   geom_point(aes(0.5, tplt_end$med[4]), 
#              size = 8, alpha = 0.025) +
#   #precip_cv_year
#   geom_jitter(data = tplt_ss_end, aes(1.5, precip_cv_year),
#               shape = 1, alpha = 0.2, size = 3) +
#   geom_segment(aes(y = tplt_end$lci[5], yend = tplt_end$uci[5],
#                    x = 1.5, xend = 1.5), 
#                lineend = "round",
#                linewidth = 3, alpha = 0.025) +
#   geom_point(aes(1.5, tplt_end$med[5]), 
#              size = 8, alpha = 0.025) +
#   xlim(-2, 2) +
#   # ylim(c(-0.28, 0.1)) +
#   ylab('temp_sd_season') +
#   xlab('') +
#   theme_bw() +
#   theme(legend.position = 'none',
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(linewidth = 2),
#         axis.ticks = element_line(size = 1.5),
#         axis.text = element_text(size = 16),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title = element_text(size = 18),
#         axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
#         axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
#         axis.ticks.length = unit(0.2, 'cm')) #length of axis tick



#effect sizes
# % increase in Gen Length for 1 unit change in sd season
(exp(summary(f1)$coefficient[3,1]) - 1) * 100
# % increase in Gen Length for 1 unit change in sd year
(exp(summary(f1)$coefficient[4,1]) - 1) * 100
# % increase in Gen Length for 1 unit change in cv precip season
(exp(summary(f1)$coefficient[5,1]) - 1) * 100
# % increase in Gen Length for 1 unit change in cv precip year
(exp(summary(f1)$coefficient[6,1]) - 1) * 100


# %increase in Gen Length over observed range of sd season
(exp(summary(f1)$coefficient[3,1] * 
       diff(range(bird_df$temp_sd_season))) - 1) * 100
# %increase in Gen Length over observed range of sd year
(exp(summary(f1)$coefficient[4,1] * 
       diff(range(bird_df$temp_sd_year))) - 1) * 100
# %increase in Gen Length over observed range of cv precip season
(exp(summary(f1)$coefficient[5,1] * 
       diff(range(bird_df$precip_cv_season))) - 1) * 100
# %increase in Gen Length over observed range of sd year
(exp(summary(f1)$coefficient[6,1] * 
       diff(range(bird_df$precip_cv_year))) - 1) * 100



#lMass
pr_fun(model = f1, 
       var_interest = 'lMass',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year',
                     'precip_cv_season',
                     'precip_cv_year'))
#temp season
pr_fun(model = f1, 
       var_interest = 'temp_sd_season',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year',
                     'precip_cv_season',
                     'precip_cv_year'))
#temp year
pr_fun(model = f1, 
       var_interest = 'temp_sd_year',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year',
                     'precip_cv_season',
                     'precip_cv_year'))
#precip season
pr_fun(model = f1, 
       var_interest = 'precip_cv_season',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year',
                     'precip_cv_season',
                     'precip_cv_year'))
#precip year
pr_fun(model = f1, 
       var_interest = 'precip_cv_year',
       var_names = c('lMass',
                     'temp_sd_season', 
                     'temp_sd_year',
                     'precip_cv_season',
                     'precip_cv_year'))


# PCA intra and inter -----------------------------------------------------

#same results as raw, essentially
tt_pca <- dplyr::select(bird_df, 
                        temp_sd_year, 
                        temp_sd_season) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

t1 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df)
t2 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df)
t3 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,2] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df)
summary(t1)
AIC(t1)
AIC(t2)
AIC(t3)

#function to backtransform from PCA 
#transform back to original data
bt_fun <- function(vals, pca_obj)
{
  t(t(vals %*% t(pca_obj$rotation)) * pca_obj$scale + pca_obj$center)
}


tt_pca$center + bt_fun(vals = cbind(1, 1), 
                       pca_obj = tt_pca)


c1 <- lme4::lmer(lGL ~ lMass + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + lMass | fac_Family) + 
                   (-1 + temp_sd_season | fac_Family) + 
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df)
summary(c1)


# PCA intra, inter, Mass -----------------------------------------------------

#same results, essentially
tt_pca <- dplyr::mutate(bird_df, 
                        lMass = lMass) %>%
  dplyr::select(lMass,
                temp_sd_year, 
                temp_sd_season) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

t1 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)
t2 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,2] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df)
t3 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,3] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)
t4 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + + 
                   (1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)
t5 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,1] +
                   (1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df)
t6 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,2] +
                   (1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df)
t7 <- lme4::lmer(lGL ~ 
                   tt_pca$x[,3] +
                   (1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)

summary(t1)
AIC(t1) #best
AIC(t2)
AIC(t3)
AIC(t4)
AIC(t5)
AIC(t6)
AIC(t7)


# PCA DHI intra inter -----------------------------------------------------

#lower var both within and across years -> longer gen time
bird_dft <- dplyr::filter(bird_df, !is.na(dhi_cv_year))
tt_pca <- dplyr::select(bird_dft,
                        dhi_cv_year, 
                        dhi_cv_season) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

t1 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_dft)
t2 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_dft)
t3 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,2] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_dft)
summary(t1)
AIC(t1) #best
AIC(t2)
AIC(t3)


c1 <- lme4::lmer(lGL ~ lMass + 
                   dhi_cv_year + 
                   dhi_cv_season + 
                   (1 + lMass | fac_Family) + 
                   (-1 + dhi_cv_year | fac_Family) + 
                   (-1 + dhi_cv_season | fac_Family), 
                 data = bird_dft)
c2 <- lme4::lmer(lGL ~ lMass + 
                   dhi_cv_year + 
                   (1 + lMass | fac_Family) + 
                   (-1 + dhi_cv_year | fac_Family), 
                 data = bird_dft)
c3 <- lme4::lmer(lGL ~ lMass + 
                   dhi_cv_season + 
                   (1 + lMass | fac_Family) + 
                   (-1 + dhi_cv_season | fac_Family), 
                 data = bird_dft)

car::vif(c1)
AIC(c1) #best
AIC(c2)
AIC(c3) #close
summary(c1)
summary(c3)


# PCA intra, inter, color -------------------------------------------------

tt_pca <- dplyr::select(bird_df, 
                        temp_sd_year, 
                        temp_sd_season, 
                        temp_sp_color_year) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

t1 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)
t2 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,1] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df)
t3 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,2] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df)
t4 <- lme4::lmer(lGL ~ lMass + 
                   tt_pca$x[,3] + 
                   (1 + lMass | fac_Family) + 
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df)
summary(t1)
AIC(t1) #best
AIC(t2)
AIC(t3)
AIC(t4)


# SEM explore -------------------------------------------------------------

semdata <- data.frame(bird_df, pc1, pc2) %>%
  dplyr::select(Mass, 
                temp_sd_season, 
                temp_sd_year,
                temp_mean,
                temp_sp_color_year,
                pc1,
                pc2,
                GenLength) %>%
  dplyr::mutate(lmass = lMass,
                lgl = lGL)

model1 <- '
  lmass ~ pc1 + pc2
  lgl ~ pc1 + pc2 + lmass'

model1.fit <- lavaan::sem(model1, data = semdata) 
summary(model1.fit)
AIC(model1.fit)





# pasted from 5-novar.R resids --------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# run_date <- '2023-09-26'
# bird_df <- readRDS(paste0(dir, '/Results/ge-bird-novar-', run_date,
#                           '/ge-bird-novar-data-', run_date, '.rds'))$pro_data
# fit <- readRDS(paste0(dir, '/Results/ge-bird-novar-', run_date,
#                       '/ge-bird-novar-fit-', run_date, '.rds'))

mu_mn <- MCMCvis::MCMCpstr(fit, params = 'mu')[[1]]
bird_df$resids <- bird_df$lGL - mu_mn

str(bird_df)
ggplot(bird_df, aes(cen_lon, cen_lat, col = resids)) +
  geom_point(size = 2) +
  theme_bw()

str(bird_df)
plot(factor(bird_df$Order), bird_df$resids)
plot(factor(bird_df$Family), bird_df$resids)
plot(bird_df$Beak.Depth, bird_df$resids)
summary(lm(bird_df$resids ~ bird_df$Beak.Depth))
plot(factor(bird_df$Order), bird_df$Beak.Depth)
plot(bird_df$Wing.Length, bird_df$resids)
summary(lm(bird_df$resids ~ bird_df$Wing.Length))

plot(factor(bird_df$Habitat), bird_df$resids)
summary(lm(bird_df$resids ~ factor(bird_df$Habitat)))
plot(factor(bird_df$Habitat.Density), bird_df$resids)
summary(lm(bird_df$resids ~ factor(bird_df$Habitat.Density)))

plot(factor(bird_df$Trophic.Level), bird_df$resids)
summary(lm(bird_df$resids ~ factor(bird_df$Trophic.Level)))
plot(factor(bird_df$Trophic.Niche), bird_df$resids)
summary(lm(bird_df$resids ~ factor(bird_df$Trophic.Niche)))
plot(factor(bird_df$Primary.Lifestyle), bird_df$resids)
plot(bird_df$temp_rng_season, bird_df$resids)
summary(lm(bird_df$resids ~ bird_df$temp_rng_season))


dplyr::group_by(bird_df, Family) %>%
  dplyr::summarize(mn_bd = mean(Beak.Depth)) %>%
  dplyr::arrange(desc(mn_bd))

dplyr::group_by(bird_df, Family) %>%
  dplyr::summarize(mn_r = mean(resids), n = n()) %>%
  dplyr::arrange(desc(mn_r)) %>%
  tail() %>%
  print(n = 50)

plot(bird_df$temp_mean, bird_df$resids)
summary(lm(bird_df$resids ~ bird_df$temp_mean))

dplyr::filter(bird_df, Habitat == 'Coastal')
dplyr::filter(bird_df, Trophic.Level == 'Scavenger')
dplyr::filter(bird_df, Trophic.Niche == 'Nectarivore')

ggplot(bird_df, aes(Wing.Length, resids, col = Order)) +
  geom_point() +
  theme_bw()

ps <- dplyr::filter(bird_df, Order == 'Psittaciformes')

ggplot(ps, aes(Beak.Depth, resids, color = Family)) +
  geom_point() +
  theme_bw()

#Brain size from Griesser et al.
bs <- read.csv(paste0(dir, 'data/L1/trait/Griesser_et_al_2023_PNAS.csv'))
tt <- dplyr::mutate(bird_df, tip_label = stringr::str_to_title(gsub(' ' , '_', Birdtree_name))) %>%
  dplyr::left_join(bs, by = 'tip_label') %>%
  dplyr::filter(!is.na(brain))

#Cooney - has clutch size but not as many as Bird et al.
co <- read.csv(paste0(dir, 'data/L1/trait/Cooney_et_al_2020_Nat_comms.csv'))
tt2 <- dplyr::mutate(bird_df, binomial = stringr::str_to_title(gsub(' ' , '_', Birdtree_name))) %>%
  dplyr::left_join(co, by = 'binomial') %>%
  dplyr::filter(!is.na(log_clutch_size))

#clutch size from Bird et al.
bb <- read.csv(paste0(dir, 'data/L1/trait/bird_et_al_clutch_size.csv'))

#survival etc from Bird et al.
bdat <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::select(-Order, -Family, -GenLength)

tt4 <- dplyr::left_join(bird_df, bdat, by = c('Birdtree_name' = 'Sci_name')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::mutate(Scientific.name = stringr::str_to_sentence(Birdtree_name),
                #correct for incorrect scaling in data
                Measured_survival = Measured_survival * 0.01) %>%
  dplyr::left_join(dplyr::select(bb, -Order, -Family, -Genus), 
                   by = 'Scientific.name') %>%
  dplyr::select(ID, Birdtree_name, Avonet_name, Family, Order, Trophic.Niche,
                temp_mean, temp_sd_year, temp_sd_season, temp_sp_color_month,
                precip_mean, precip_cv_year, precip_cv_season, precip_sp_color_month,
                Mass, GenLength, lMass, lGL, Mean.clutch.size, Measured_survival, 
                Measured_age_first_breeding, Measured_max_longevity, 
                Modeled_survival, Modeled_age_first_breeding, Modeled_max_longevity)

#PCA
#105 species
pdat <- dplyr::filter(tt4, 
                      #!is.na(Mean.clutch.size),
                      !is.na(Measured_survival),
                      !is.na(Measured_age_first_breeding),
                      !is.na(Measured_max_longevity))
NROW(pdat)

#don't know why survival values are so high for measured (this is how it is in paper)
plot(tt4$Measured_survival, tt4$Modeled_survival)
abline(0, 1, lty = 2, lwd = 2)
plot(tt4$Measured_age_first_breeding, tt4$Modeled_age_first_breeding)
plot(tt4$Measured_max_longevity, tt4$Modeled_max_longevity)

#3011 species
pdat <- dplyr::filter(tt4, 
                      #!is.na(Mean.clutch.size),
                      !is.na(Modeled_survival),
                      !is.na(Modeled_age_first_breeding),
                      !is.na(Modeled_max_longevity))
NROW(pdat)

tt_pca <- pdat[,c(#'Mean.clutch.size', 
  'Measured_survival',
  'Measured_age_first_breeding',
  'Measured_max_longevity')] %>% 
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

tt_pca <- pdat[,c(#'Mean.clutch.size', 
                  'Modeled_survival',
                  'Modeled_age_first_breeding',
                  'Modeled_max_longevity')] %>% 
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

#essentially the same thing as GL
plot(pdat$lGL, tt_pca$x[,1])

f1 <- lm(lGL ~ #lMass + 
           temp_sd_season + 
           temp_sd_year +
           temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year +
         precip_sp_color_month,
         #brain,
         data = pdat)
summary(f1)
car::vif(f1)

f1b <- lm(lGL ~ #lMass + 
           temp_sd_season + 
           temp_sd_year +
           temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year +
           precip_sp_color_month,
         #brain,
         data = tt4)
summary(f1b)
car::vif(f1b)

#only species with data (no total interps)
pd3 <- dplyr::filter(tt4, !is.na(Measured_survival) |
                       !is.na(Measured_age_first_breeding) |
                       !is.na(Measured_max_longevity))
#r2 = 0.013
f1c1 <- lm(lGL ~ #lMass + 
            temp_sd_season + 
            temp_sd_year +
            temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year +
            precip_sp_color_month,
          #brain,
          data = pd3)
summary(f1c1)
car::vif(f1c1)

#r2 = 0.2
f2 <- lm(log(Mean.clutch.size) ~ lMass + 
           temp_sd_season + 
           temp_sd_year +
           # temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year,# +
         # precip_sp_color_month,
         #brain,
         data = tt4)
summary(f2)
car::vif(f2)

#r2 = 0.06, 0.04
f3 <- lm(Measured_survival ~ lMass +
# f3 <- lm(Modeled_survival ~ #lMass +
           temp_sd_season + 
           temp_sd_year +
           #temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year, #+
           #precip_sp_color_month,
         #brain,
         data = idm2)
summary(f3)
car::vif(f3)

f4 <- lm(Measured_age_first_breeding ~ #lMass + 
# f4 <- lm(Modeled_age_first_breeding ~ lMass + 
           temp_sd_season + 
           temp_sd_year +
           #temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year, #+
         #precip_sp_color_month,
         #brain,
         data = tt4)
summary(f4)
car::vif(f4)

#r2 = 0.03
f5 <- lm(Measured_max_longevity ~ #lMass + 
# f5 <- lm(Modeled_max_longevity ~ lMass + 
           temp_sd_season + 
           temp_sd_year +
           #temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year, #+
         #precip_sp_color_month,
         #brain,
         data = tt4)
summary(f5)
car::vif(f5)



str(bird_df)
ggplot(bird_df, aes(temp_sd_season, lGL, color = factor(Trophic.Niche))) +
  # ggplot(bird_df, aes(temp_sd_year, lGL, color = factor(Trophic.Niche))) +
  # ggplot(bird_df, aes(temp_sp_color_month, lGL, color = factor(Trophic.Niche))) +
  # ggplot(bird_df, aes(precip_cv_season, lGL, color = factor(Trophic.Niche))) +
  # ggplot(bird_df, aes(precip_cv_year, lGL, color = factor(Trophic.Niche))) +
  # ggplot(bird_df, aes(precip_sp_color_month, lGL, color = factor(Trophic.Niche))) +
  geom_point(alpha = 0.5) + 
  geom_line(stat = 'smooth', alpha = 0.5, linewidth = 2,
            method = 'lm') +
  theme_bw()


str(tt4)

# IMPUTE TRAITS -----------------------------------------------------------

library(tidyverse)
library(ape)
library(Rphylopars)
library(rstan)
library(MCMCvis)

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
run_date <- '2023-09-26'
bird_df <- readRDS(paste0(dir, '/Results/ge-bird-novar-', run_date,
                          '/ge-bird-novar-data-', run_date, '.rds'))$pro_data

#clutch size from Bird et al.
bb <- read.csv(paste0(dir, 'data/L1/trait/bird_et_al_clutch_size.csv'))

#survival etc from Bird et al.
bdat <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::select(-Order, -Family, -GenLength)

tt4 <- dplyr::left_join(bird_df, bdat, by = c('Birdtree_name' = 'Sci_name')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::mutate(Scientific.name = stringr::str_to_sentence(Birdtree_name),
                #correct for incorrect scaling in data
                Measured_survival = Measured_survival * 0.01) %>%
  dplyr::left_join(dplyr::select(bb, -Order, -Family, -Genus), 
                   by = 'Scientific.name') %>%
  dplyr::select(ID, Birdtree_name, Avonet_name, Family, Order, Trophic.Niche,
                temp_mean, temp_sd_year, temp_sd_season, temp_sp_color_month,
                precip_mean, precip_cv_year, precip_cv_season, precip_sp_color_month,
                Mass, GenLength, lMass, lGL, Mean.clutch.size, Measured_survival, 
                Measured_age_first_breeding, Measured_max_longevity, 
                Modeled_survival, Modeled_age_first_breeding, Modeled_max_longevity)

dd <- dplyr::mutate(tt4,
              species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::mutate(Measured_log_age_first_breeding = log(Measured_age_first_breeding),
                Measured_log_max_longevity = log(Measured_max_longevity),
                Measured_log_clutch_size = log(Mean.clutch.size)) %>%
  dplyr::select(species, 
                lMass,
                Measured_survival,
                Measured_log_age_first_breeding,
                Measured_log_max_longevity,
                Measured_log_clutch_size)

#read in phylo tree (birdtree.org - Ericson 0001-1000)
tr <- ape::read.tree(paste0(dir, 'data/L1/trait/AllBirdsEricson1.tre'))

#df with names and idx
idx_df <- data.frame(idx = 1:NROW(dd), 
                     name = stringr::str_to_title(gsub(' ', '_', dd$species)))

#species not found in both datasets (species to drop from tree)
nm <- setdiff(tr[[1]]$tip.label, idx_df$name)

#prune specified tips from all trees
pr_tr <- lapply(tr, ape::drop.tip, tip = nm)
class(pr_tr) <- "multiPhylo"
tree_n <- pr_tr[[1]]

#get index for name order on tips
j_idx <- dplyr::left_join(data.frame(name = tree_n$tip.label), idx_df, 
                          by = 'name')

#run phylo imputation
ir <- Rphylopars::phylopars(trait_data = dd[j_idx$idx,], 
                            tree = tree_n, 
                            phylo_correlated = TRUE,
                            # model = "BM") #20804
                            # model = "OU") #1942203
                            model = "lambda") #10755
                            # model = "mvOU") # threw an error
                            # model = "delta") # threw an error
                            # model = "EB") #20697
                            # model = "star")

AIC(ir)

#returns species means as well as some internal nodes
#get just species rows
ridx <- which(row.names(ir$anc_rec) %in% dd$species)

#variance (uncertainty)
iv <- data.frame(species = dd$species[j_idx$idx], 
                 ir$anc_var[ridx,]) %>%
  dplyr::mutate(SD_survival = sqrt(Measured_survival),
                SD_log_age_first_breeding = sqrt(Measured_log_age_first_breeding),
                SD_log_max_longevity = sqrt(Measured_log_max_longevity),
                SD_log_clutch_size = sqrt(Measured_log_clutch_size)) %>%
  dplyr::arrange(species) %>%
  dplyr::select(species, SD_survival, 
                SD_log_age_first_breeding,
                SD_log_max_longevity,
                SD_log_clutch_size)


id2 <- data.frame(species = dd$species[j_idx$idx], 
                    ir$anc_rec[ridx,]) %>%
  dplyr::arrange(species) %>%
  dplyr::rename(Phylo_survival = Measured_survival,
                Phylo_log_age_first_breeding = Measured_log_age_first_breeding,
                Phylo_log_max_longevity = Measured_log_max_longevity,
                Phylo_log_clutch_size = Measured_log_clutch_size) %>%
  #join with measured ta
  dplyr::left_join(dplyr::select(dd, -lMass), by = 'species') %>%
  #join with var
  dplyr::left_join(iv, by = 'species')
row.names(id) <- NULL

#run 
# Measured_survival ~ mass + env
# Measured_age_first_breeding ~ mass + env
# Measured_max_longevity ~ mass + env
# Mean.clutch.size ~ mass + env

#join with env data
idm2 <- dplyr::mutate(tt4,
              species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::select(species,
                temp_mean,
                temp_sd_year,
                temp_sd_season,
                temp_sp_color_month,
                precip_cv_year,
                precip_cv_season,
                precip_sp_color_month,
                Modeled_survival,
                Modeled_age_first_breeding,
                Modeled_max_longevity) %>%
  dplyr::left_join(id2, by = 'species')


#models

#0.17, 0.01 -> lambda
#0.09, 0.02
ti1 <- lm(Phylo_survival ~ lMass + 
           temp_sd_season + 
           temp_sd_year +
           # temp_sp_color_month +
           precip_cv_season + 
           precip_cv_year,# + 
            # precip_sp_color_month,
         data = idm2)
summary(ti1)
car::vif(ti1)

ti2 <- lm(Phylo_log_age_first_breeding ~ lMass + 
            temp_sd_season + 
            temp_sd_year +
            temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year + 
            precip_sp_color_month,
          data = idm2)
summary(ti2)
car::vif(ti2)

ti3 <- lm(Phylo_log_max_longevity ~ #lMass + 
            temp_sd_season + 
            temp_sd_year +
            temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year + 
            precip_sp_color_month,
          data = idm2)
summary(ti3)
car::vif(ti3)

#0.14, 0.14
#0.14, 0.14
ti4 <- lm(Phylo_log_clutch_size ~ lMass + 
            temp_sd_season + 
            temp_sd_year +
            # temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year, #+ 
            # precip_sp_color_month,
          data = idm2)
summary(ti4)
car::vif(ti4)

#test resids for phylo signal
idm2 %>% print(width = Inf)
#apply to residuals
resid_srt <- residuals(ti1)[j_idx$idx]
resid_srt <- residuals(ti4)[j_idx$idx]
names(resid_srt) <- j_idx$name
phy_res <- picante::phylosignal(resid_srt, tree_n) #quite slow

########
#subset to test model speed of PGLS
# set.seed(1)
# sidx <- sample(1:NROW(idm2), size = 3000)
# idm_s <- idm2[sidx,]
idm_s <- idm2
#prune specified tips from all trees
nm2 <- setdiff(tr[[1]]$tip.label, idm_s$species)
pr_tr2 <- lapply(tr, ape::drop.tip, tip = nm2)
class(pr_tr2) <- "multiPhylo"
tree_n2 <- pr_tr2[[1]]
#df with names and idx
idx_df2 <- data.frame(idx = 1:NROW(idm_s), 
                     name = stringr::str_to_title(gsub(' ', '_', idm_s$species)))
#get index for name order on tips
j_idx2 <- dplyr::left_join(data.frame(name = tree_n2$tip.label), idx_df2, 
                          by = 'name')
idm_s2 <- idm_s[j_idx2$idx,]
#run phylo imputation
library(nlme)
pglsModel <- gls(Phylo_survival ~ lMass +
                   temp_sd_season +
                   temp_sd_year +
                   precip_cv_season +
                   precip_cv_year,
                 correlation = ape::corBrownian(phy = tree_n2),
                 data = idm_s2,
                 method = "ML")
summary(pglsModel)

library(caper)
comp.data <- caper::comparative.data(tree_n2, 
                                     as.data.frame(idm_s2), 
                                     names.col = 'species',
                                     vcv = TRUE)
pf1 <- caper::pgls(Phylo_survival ~ lMass + 
                  temp_sd_season + 
                  temp_sd_year +
                  precip_cv_season + 
                  precip_cv_year, 
                data = comp.data,
                lambda = 'ML')
summary(pf1)

pf2 <- caper::pgls(Phylo_log_clutch_size ~ lMass + 
                     temp_sd_season + 
                     temp_sd_year +
                     precip_cv_season + 
                     precip_cv_year, 
                   data = comp.data,
                   lambda = 'ML')
summary(pf2)
########

hist(idm2$Phylo_survival)
hist(idm2$SD_survival)
hist(idm2$Phylo_log_clutch_size)
hist(idm2$SD_log_clutch_size)

#pca for LH traits
tt_pca <- idm2[,c('Phylo_log_clutch_size', 
                 'Phylo_survival')] %>%
                 # 'lMass')] %>%
  # 'Phylo_survival',
  # 'Phylo_log_age_first_breeding',
  # 'Phylo_log_max_longevity')] %>%
  prcomp(center = TRUE, scale. = TRUE)

factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

#THIS IS PROBS THE MOST USEFUL
factoextra::fviz_pca_biplot(tt_pca,
                            axes = c(1,2), 
                            label = "var",
                             habillage = tt4$Order,
                             addEllipses = TRUE, 
                             ellipse.level = 0.95)

idm2[which(tt4$Order == 'Anseriformes'),] %>%
  dplyr::select(species, Measured_survival, Measured_age_first_breeding, 
                Measured_max_longevity, Mean.clutch.size)
hist(idm$Measured_survival)
hist(idm$Mean.clutch.size)

#0.14, 0.08 -> surv vs. egg
ti5 <- lm(tt_pca$x[,1] ~ lMass + 
            temp_sd_season + 
            temp_sd_year +
            # temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year,# + 
          # precip_sp_color_month,
          data = idm2)
summary(ti5)
car::vif(ti5)
#less mass = higher PC = more eggs/lower surv
#more sd season = higher PC = more eggs/lower surv  -> supports what we know
#more sd year = lower PC = fewer eggs/higher surv   -> supports expectations 
AIC(ti5)

#0.52, 0.05
#0.17, 0.06 -> both + vs -
ti6 <- lm(tt_pca$x[,2] ~ lMass + 
            temp_sd_season + 
            temp_sd_year +
            # temp_sp_color_month +
            precip_cv_season + 
            precip_cv_year, #+ 
          # precip_sp_color_month,
          data = idm2)
summary(ti6)
car::vif(ti6)
#more mass = higher PC = more eggs/higher surv
#more sd season = higher PC = more eggs/higher surv  -> supports expectations
#more sd year = higher PC = more eggs/higher surv    -> supports expectations

