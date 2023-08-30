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
library(viridis)


# read in data -------------------------------------------------

bird_df <- read.csv(paste0(dir, 'Data/L2/main-bird-data.csv')) %>%
  dplyr::filter(!is.na(GenLength)) %>%
  dplyr::mutate(fac_Family = factor(Family),
                fac_Order = factor(Order))

#filter - no marine/coastal, only residents
'%ni%' <- Negate('%in%')
bird_df2 <- dplyr::filter(bird_df, 
                          Habitat %ni% c('Marine', 'Coastal'),
                          Primary.Lifestyle %ni% ('Aquatic'),
                          Migration == 1)


# plot gen length on map ----------------------------------------------

#should replace this with rasterized gen length on map

#map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#base plot
base_plt <- ggplot() + 
  #plot map
  geom_sf(data = world,
          #fill = 'grey89') +
          fill = 'aliceblue') +
  # theme_bw()
  #graticles
  theme(panel.grid.major = #element_line(color = black,
                            #            linetype = 'dashed', 
                             #           linewidth = 0.1),
          element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = 'NA')) +
  scale_x_continuous(breaks = seq(-180, 180, by = 15))+
  scale_y_continuous(breaks = seq(-80, 80, by = 20))

#pts
gen_pts_wo_seabirds_coastal <- base_plt +
  #plot pts
  geom_point(data = bird_df2,
             aes(cen_lon, cen_lat,
                 col = log(GenLength)),
             alpha = 0.5,
             size = 1.1) +
  scale_color_viridis()

gen_pts_all <- base_plt +
  #plot pts
  geom_point(data = bird_df,
             aes(cen_lon, cen_lat,
                 col = log(GenLength)),
             alpha = 0.5,
             size = 1.1) +
  scale_color_viridis()

gen_pts_wo_seabirds_coastal
gen_pts_all

#plot by habitat type
(uhab <- unique(bird_df$Habitat))
hab_fun <- function(habitat)
{
  bird_t <- dplyr::filter(bird_df, Habitat == habitat)
  fun_plt <- base_plt +
    geom_point(data = bird_t,
               aes(cen_lon, cen_lat,
                   col = log(GenLength)),
               alpha = 0.8,
               size = 1.1) +
    scale_color_viridis() +
    ggtitle(habitat) +
    xlab('Longitude') +
    ylab('Latitude')

    print(fun_plt)
}

hab_fun(habitat = 'Forest')
hab_fun(habitat = 'Woodland')
hab_fun(habitat = 'Wetland')
hab_fun(habitat = 'Shrubland')
hab_fun(habitat = 'Grassland')
hab_fun(habitat = 'Rock')
hab_fun(habitat = 'Desert')
hab_fun(habitat = 'Coastal')
hab_fun(habitat = 'Marine')
hab_fun(habitat = 'Riverine')
hab_fun(habitat = 'Human Modified')


# quick and dirty GAM spatial model ---------------------------------------

# #variation over space
# sg_fit <- mgcv::gam(log(bird_df2$GenLength) ~ s(bird_df2$cen_lon, 
#                                                            bird_df2$cen_lat))
# 
# plot(sg_fit, contour.col = 'black', too.far = 0.2, scheme = 2, rug = TRUE,
#      main = 'GAM spatial - Red = short gen, Yellow = long gen',
#      xlab = 'Longitude',
#      ylab = 'Latitude')


# mass ~ env --------------------------------------------------------------

f1 <- lme4::lmer(log(Mass) ~ env1_pc1 + 
                   env1_pc2 + 
                   env1_pc3 + 
                   (-1 + env1_pc1 | fac_Family) + 
                   (-1 + env1_pc2 | fac_Family) +
                   (-1 + env1_pc3 | fac_Family), 
                 data = bird_df2)

summary(f1)


# how long before temp trend exceeds 2 sd interannual ---------------------

plot(bird_df2$temp_sd_year, bird_df2$temp_slope, 
     col = rgb(0,0,0,0.1),
     pch = 19)

#ex = number of years before temp exceeds 2 sd
#AND
#delta_t = how much temp will change in 1 generation (in sds)
#AND
#delta_haldane = how much change (in sd) per generation
#AND
#n_gen = how many gens before temp will exceed 2 sd
ex_df <- dplyr::mutate(bird_df2, 
                       ex = 2 * temp_sd_year / temp_slope,
                       delta_t = temp_slope / temp_sd_year * GenLength,
                       #(degrees / year) * (sd / degrees) * (year / gen) = sd / gen
                       delta_haldane = (temp_slope / temp_sd_year) * GenLength,
                       n_gen = ex / GenLength)
hist(ex_df$ex, xlim = c(0, 125), breaks = 10000)
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

#Bürger and Lynch 1995: considerably less than 10% per gen
#Gingerich 2009: 0.1 to 0.3 SD per generation


# visualize ---------------------------------------------------------------

# plt <- dplyr::select(tt_env3, lon, lat, env3_pc1) %>%
#   # dplyr::mutate(tt = scale(temp_sd_resid, scale = TRUE)) %>%
#   # dplyr::mutate(tt = log(temp_sd_resid)) %>%
#   dplyr::mutate(tt = env3_pc1) %>%
#   dplyr::select(lon, lat, tt) %>%
#   terra::rast(crs = "epsg:4326")
# rv <- range(terra::values(plt$tt), na.rm = TRUE)
# pal <- leaflet::colorNumeric(palette = "RdBu",
#                              domain = rv,
#                              reverse = T)
# plot(plt, range = rv, col = pal(seq(rv[1], rv[2], by = 0.1)))


# raw temp -----------------------------------------------------------

#candidate models
f1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sp_color_year + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sp_color_year | fac_Family) + 
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
f2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sp_color_year + 
                   temp_sd_season + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sp_color_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family), 
                 data = bird_df2)
f3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sp_color_year + 
                   temp_sd_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sp_color_year | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
f4 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_season + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_season | fac_Family), 
                 data = bird_df2)
f5 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
f6 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
f7 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   (1 + log(Mass) | fac_Family), 
                 data = bird_df2)
summary(f1)
summary(f4)
summary(f5)
summary(f6) #best (sd season and sd year and Mass)

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
# % increase in Gen Length for 1 unit change in sd season
(exp(summary(f6)$coefficient[3,1]) - 1) * 100
# % increase in Gen Length for 1 unit change in sd year
(exp(summary(f6)$coefficient[4,1]) - 1) * 100
# %increase in Gen Length over observed range of sd season
(exp(summary(f6)$coefficient[3,1] * 
       diff(range(bird_df2$temp_sd_season))) - 1) * 100
# %increase in Gen Length over observed range of sd year
(exp(summary(f6)$coefficient[4,1] * 
       diff(range(bird_df2$temp_sd_year))) - 1) * 100

#model selection
AIC(f1)
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6) #best (no color)
AIC(f7)
MuMIn::r.squaredGLMM(f1)
MuMIn::r.squaredGLMM(f2)
MuMIn::r.squaredGLMM(f3)
MuMIn::r.squaredGLMM(f5)


#VIF - looks good
car::vif(f1)
vif_1 <- lme4::lmer(log(Mass) ~
                      temp_sd_season + 
                      temp_sd_year + 
                      (1 + temp_sd_season | fac_Family) +
                      (-1 + temp_sd_year | fac_Family), 
                    data = bird_df2)
vif_2 <- lme4::lmer(temp_sd_season ~
                      log(Mass) + 
                      temp_sd_year + 
                      (1 + log(Mass) | fac_Family) +
                      (-1 + temp_sd_year | fac_Family), 
                    data = bird_df2)
vif_3 <- lme4::lmer(temp_sd_year ~
                      log(Mass) +
                      temp_sd_season + 
                      (1 + log(Mass) | fac_Family) +
                      (-1 + temp_sd_season | fac_Family), 
                    data = bird_df2)

vif_1 <- summary(lm(log(Mass) ~ temp_sd_year + temp_sd_season, data = bird_df2))
vif_2 <- summary(lm(temp_sd_year ~ log(Mass) + temp_sd_season, data = bird_df2))
vif_3 <- summary(lm(temp_sd_season ~ log(Mass) + temp_sd_year, data = bird_df2))
vif_4 <- summary(lm(temp_sd_year ~ temp_sd_season, data = bird_df2))
1 / (1 - vif_1$adj.r.squared)
1 / (1 - vif_2$adj.r.squared)
1 / (1 - vif_3$adj.r.squared)
1 / (1 - vif_4$adj.r.squared)
1 / (1 - MuMIn::r.squaredGLMM(vif_1)[,1])
1 / (1 - MuMIn::r.squaredGLMM(vif_2)[,1])
1 / (1 - MuMIn::r.squaredGLMM(vif_3)[,1])


#partial residual plots
library(remef)
#effect of Mass
coef_fit <- coef(summary(f6))
yp <- remef::remef(f6, fix = c('temp_sd_season', 
                               'temp_sd_year'),
                   ran = "all")
plot(log(bird_df2$Mass), yp, col = rgb(0,0,0,0.1), pch = 19)
abline(a = coef_fit[1,1], b = coef_fit[2,1], col = 'red')
#effect of seasonality
yp <- remef::remef(f6, fix = c('log(Mass)', 
                               'temp_sd_year'),
                   ran = "all")
plot(bird_df2$temp_sd_season, yp, col = rgb(0,0,0,0.1), pch = 19)
abline(a = coef_fit[1,1], b = coef_fit[3,1], col = 'red')
#effect of inter-annual var
yp <- remef::remef(f6, fix = c('log(Mass)', 
                               'temp_sd_season'),
                   ran = "all")
plot(bird_df2$temp_sd_year, yp, col = rgb(0,0,0,0.1), pch = 19)
abline(a = coef_fit[1,1], b = coef_fit[4,1], col = 'red')

yp <- remef::remef(f6, fix = c('log(Mass)'),
                   ran = "all")

bird_df2$yp <- yp
ggplot(bird_df2, aes(temp_sd_year, temp_sd_season, col = yp)) +
  geom_point()
summary(lm(yp ~ temp_sd_year + temp_sd_season, data = bird_df2))




# raw temp and precip -----------------------------------------------------

#candidate models
f1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                           temp_sd_season + 
                           temp_sd_year +
                           precip_cv_season + 
                           precip_cv_year +
                           (1 + log(Mass) | fac_Family) +
                           (-1 + temp_sd_season | fac_Family) +
                           (-1 + temp_sd_year | fac_Family) +
                           (-1 + precip_cv_season | fac_Family) +
                           (-1 + precip_cv_year | fac_Family), 
                         data = bird_df2)
f2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   precip_cv_season + 
                   precip_cv_year +
                   (1 + log(Mass) | fac_Family) +
                   (-1 + precip_cv_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df2)
f3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_season + 
                   temp_sd_year +
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
f4 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = bird_df2)
f5 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_year + 
                   precip_cv_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df2)
f6 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_year + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_year | fac_Family), 
                 data = bird_df2)
f7 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_year + 
                   temp_sd_season + 
                   precip_cv_season + 
                   (1 + log(Mass) | fac_Family) +
                   (-1 + temp_sd_year | fac_Family) +
                   (-1 + temp_sd_season | fac_Family) +
                   (-1 + precip_cv_season | fac_Family), 
                 data = bird_df2)
summary(f1)

AIC(f1) #best
AIC(f2)
AIC(f3)
AIC(f4)
AIC(f5)
AIC(f6)
AIC(f7)

summary(f1)
summary(f4)
summary(f5)
summary(f6) #best (sd season and sd year and Mass)



# PCA intra and inter -----------------------------------------------------

#same results as raw, essentially
tt_pca <- dplyr::select(bird_df2, 
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

t1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df2)
t2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df2)
t3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,2] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df2)
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


c1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   temp_sd_season + 
                   temp_sd_year + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + temp_sd_season | fac_Family) + 
                   (-1 + temp_sd_year | fac_Family), 
                 data = bird_df2)
summary(c1)


# PCA intra, inter, Mass -----------------------------------------------------

#same results, essentially
tt_pca <- dplyr::mutate(bird_df2, 
                        lMass = log(Mass)) %>%
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

t1 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)
t2 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,2] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df2)
t3 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,1] +
                   tt_pca$x[,3] + 
                   (1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)
t4 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + + 
                   (1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)
t5 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,1] +
                   (1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df2)
t6 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,2] +
                   (1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df2)
t7 <- lme4::lmer(log(GenLength) ~ 
                   tt_pca$x[,3] +
                   (1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)

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
bird_dft <- dplyr::filter(bird_df2, !is.na(dhi_cv_year))
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

t1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_dft)
t2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_dft)
t3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,2] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_dft)
summary(t1)
AIC(t1) #best
AIC(t2)
AIC(t3)


c1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   dhi_cv_year + 
                   dhi_cv_season + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + dhi_cv_year | fac_Family) + 
                   (-1 + dhi_cv_season | fac_Family), 
                 data = bird_dft)
c2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   dhi_cv_year + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + dhi_cv_year | fac_Family), 
                 data = bird_dft)
c3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   dhi_cv_season + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + dhi_cv_season | fac_Family), 
                 data = bird_dft)

car::vif(c1)
AIC(c1) #best
AIC(c2)
AIC(c3) #close
summary(c1)
summary(c3)


# PCA intra, inter, color -------------------------------------------------

tt_pca <- dplyr::select(bird_df2, 
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

t1 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   tt_pca$x[,2] + 
                   tt_pca$x[,3] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family) +
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)
t2 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,1] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,1] | fac_Family), 
                 data = bird_df2)
t3 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,2] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,2] | fac_Family), 
                 data = bird_df2)
t4 <- lme4::lmer(log(GenLength) ~ log(Mass) + 
                   tt_pca$x[,3] + 
                   (1 + log(Mass) | fac_Family) + 
                   (-1 + tt_pca$x[,3] | fac_Family), 
                 data = bird_df2)
summary(t1)
AIC(t1) #best
AIC(t2)
AIC(t3)
AIC(t4)


# SEM explore -------------------------------------------------------------

semdata <- data.frame(bird_df2, pc1, pc2) %>%
  dplyr::select(Mass, 
                temp_sd_season, 
                temp_sd_year,
                temp_mean,
                temp_sp_color_year,
                pc1,
                pc2,
                GenLength) %>%
  dplyr::mutate(lmass = log(Mass),
                lgl = log(GenLength))

model1 <- '
  lmass ~ pc1 + pc2
  lgl ~ pc1 + pc2 + lmass'

model1.fit <- lavaan::sem(model1, data = semdata) 
summary(model1.fit)
AIC(model1.fit)

