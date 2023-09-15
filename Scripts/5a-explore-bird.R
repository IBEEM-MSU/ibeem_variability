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
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data.csv')) %>%
  dplyr::arrange(Accepted_name) %>%
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
plot(log(bird_ras3[[c('median_gl', 'median_dh')]]))
plot(bird_ras[['n_sp']])


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
#model selection
AIC(f1) #tie best (color)
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6) #tie best (no color)
AIC(f7)
MuMIn::r.squaredGLMM(f6)

summary(f1)
summary(f4)
summary(f5)
summary(f6)

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


#VIF - looks good (VIF < 5 is OK)
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
                   #temp_sp_color_month +
                   precip_cv_season + 
                   precip_cv_year +
                   #precip_sp_color_month +
                   (1 + lMass | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   #(-1 + temp_sp_color_month | fac_Family) +
                   (-1 + precip_cv_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), #+
                 #(-1 + precip_sp_color_month | fac_Family),
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
round(summary(f1)$coefficients, 2)

AIC(f1) #best (sp color temp is within 3)
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6)
AIC(f7)

summary(f1)
summary(f4)
summary(f5)
summary(f6)


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

