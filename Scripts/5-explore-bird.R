################
# Initial exploration bird
#
# Notes:
# - might be good to do PCA to collapse env var metrics -> e.g., mean temp, sd temp, rho_l1_temp (see end of script)
# - responses are likely to vary by habitat, trophic level, migratory status
# - might be good to get number of species coexists at a given cell? (via rasterizing range maps)
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data -------------------------------------------------

bird_df <- read.csv(paste0(dir, 'Data/L2/main-bird-data.csv'))
#add cv precip
bird_df$cv_precip <- bird_df$sd_resid_precip / bird_df$mean_precip


# plot gen grouped by habitat, trophic level, mig status -----------------------

plot(factor(bird_df$avonet.Habitat), log(bird_df$bird.et.al.GenLength),
     xlab = 'Habitat type', ylab = 'log(Gen Length)')

plot(factor(bird_df$avonet.Trophic.Level), log(bird_df$bird.et.al.GenLength),
     xlab = 'Trophic level', ylab = 'log(Gen Length)')

plot(factor(bird_df$avonet.Migration), log(bird_df$bird.et.al.GenLength),
     xlab = 'Migratory status', ylab = 'log(Gen Length)')


# plot gen length on map ----------------------------------------------

#filter - no marine/coastal
'%ni%' <- Negate('%in%')
bird_df2 <- dplyr::filter(bird_df, avonet.Habitat %ni% c('Marine', 'Coastal'))

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
  theme(panel.grid.major = element_line(color = gray(0.8),
                                        linetype = 'dashed', 
                                        linewidth = 0.5),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = 'NA')) +
  scale_x_continuous(breaks = seq(-180, 180, by = 15))+
  scale_y_continuous(breaks = seq(-80, 80, by = 20))

#pts
gen_pts_wo_seabirds_coastal <- base_plt +
  #plot pts
  geom_point(data = bird_df2,
             aes(avonet.Centroid.Longitude, avonet.Centroid.Latitude,
                 col = log(bird.et.al.GenLength)),
             alpha = 0.8,
             size = 1.1) +
  scale_color_viridis()

gen_pts_all <- base_plt +
  #plot pts
  geom_point(data = bird_df,
             aes(avonet.Centroid.Longitude, avonet.Centroid.Latitude,
                 col = log(bird.et.al.GenLength)),
             alpha = 0.8,
             size = 1.1) +
  scale_color_viridis()

gen_pts_wo_seabirds_coastal
gen_pts_all

#plot by habitat type
(uhab <- unique(bird_df$avonet.Habitat))
hab_fun <- function(habitat)
{
  bird_t <- dplyr::filter(bird_df, avonet.Habitat == habitat)
  fun_plt <- base_plt +
    geom_point(data = bird_t,
               aes(avonet.Centroid.Longitude, avonet.Centroid.Latitude,
                   col = log(bird.et.al.GenLength)),
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

dplyr::filter(bird_df, 
              avonet.Habitat == 'Grassland', 
              avonet.Centroid.Latitude > 60,
              avonet.Centroid.Longitude < 20,
              avonet.Centroid.Longitude > -20)


# quick and dirty GAM spatial model ---------------------------------------

#variation over space
sg_fit <- mgcv::gam(log(bird_df2$bird.et.al.GenLength) ~ s(bird_df2$avonet.Centroid.Longitude, 
                                                           bird_df2$avonet.Centroid.Latitude))

plot(sg_fit, contour.col = 'black', too.far = 0.2, scheme = 2, rug = TRUE,
     main = 'GAM spatial - Red = short gen, Yellow = long gen',
     xlab = 'Longitude',
     ylab = 'Latitude')

#with sd temp included as covariate
sg_fit2 <- mgcv::gam(log(bird_df2$bird.et.al.GenLength) ~ s(bird_df2$avonet.Centroid.Longitude, 
                                                           bird_df2$avonet.Centroid.Latitude) +
                       bird_df2$sd_resid_temp)

plot(sg_fit2, contour.col = 'black', too.far = 0.2, scheme = 2, rug = TRUE,
     main = 'GAM spatial + sd temp - Red = short gen, Yellow = long gen',
     xlab = 'Longitude',
     ylab = 'Latitude')
summary(sg_fit2)


# filter data -----------------------------------------------------------------

str(bird_df)

#number of species each grouping
table(bird_df$avonet.Trophic.Level)
#Migration (1 - sedentary, 2 - partially migratory, 3 - migratory)
table(bird_df$avonet.Migration)
table(bird_df$avonet.Habitat)


#filter species
tt <- dplyr::filter(bird_df, 
                    avonet.Migration == 1, #sedentary
                    # avonet.Migration == 2, #partially migratory
                    # avonet.Migration == 3, #migratory
                    avonet.Trophic.Level == 'Herbivore',
                    !is.na(bird.et.al.GenLength),
                    !is.na(sd_resid_temp))


# explore gen length ~ sd temp --------------------------------------------------

#by habitat type
ggplot(tt, aes(sd_resid_temp, log(bird.et.al.GenLength),
                    col = factor(avonet.Habitat))) +
  geom_point(alpha = 0.1) +
  geom_line(stat = 'smooth',
            method = 'lm',
            linewidth = 1.5,
            alpha = 0.4) +
  ylim(c(0, 3)) +
  theme_bw()


#separate lm for each habitat type
uhab <- unique(bird_df$avonet.Habitat)

par(mfrow = c(3,3))
out <- data.frame(habitat = rep(NA, length(uhab)),
                  N = NA,
                  sl = NA,
                  pval = NA,
                  r2 = NA)
for (i in 1:length(uhab))
{
  #i <- 1
  tdf <- dplyr::filter(tt, avonet.Habitat == uhab[i])

  out$habitat[i] <- uhab[i]
  out$N[i] <- NROW(tdf)
  
  if (NROW(tdf) > 1)
  {
    t1 <- lm(log(tdf$bird.et.al.GenLength) ~ tdf$sd_resid_temp)
    plot(tdf$sd_resid_temp, log(tdf$bird.et.al.GenLength),
         main = paste0(uhab[i]),
         pch = 19, col = rgb(0,0,0,0.1))
    
    # t1 <- lm(log(tdf$bird.et.al.GenLength) ~ tdf$cv_precip)
    # plot(tdf$cv_precip, log(tdf$bird.et.al.GenLength),
    #      main = paste0(uhab[i]),
    #      pch = 19, col = rgb(0,0,0,0.1))
    
    #summary model fit
    st1 <- summary(t1)
    
    #fill df
    out$sl[i] <- round(st1$coefficients[2,1], 3)
    out$pval[i] <- round(st1$coefficients[2,4], 3)
    out$r2[i] <- round(st1$r.squared, 3)
  }
}

out


# positive examples -------------------------------------------------------

tt2 <- dplyr::filter(bird_df2, 
                     avonet.Habitat == 'Human Modified',
                     !is.na(bird.et.al.GenLength),
                     !is.na(sd_resid_temp))

ggplot(tt2, aes(sd_resid_temp, log(bird.et.al.GenLength),
               col = factor(avonet.Trophic.Level))) +
  geom_point(alpha = 0.5) +
  geom_line(stat = 'smooth',
            method = 'lm',
            linewidth = 1.5,
            alpha = 0.8) +
  ylim(c(0, 3)) +
  theme_bw() +
  ggtitle('Human Modified birds')

tt3 <- dplyr::filter(bird_df2, 
                     avonet.Habitat == 'Grassland',
                     !is.na(bird.et.al.GenLength),
                     !is.na(sd_resid_temp))

ggplot(tt3, aes(sd_resid_temp, log(bird.et.al.GenLength),
                col = factor(avonet.Trophic.Level))) +
  geom_point(alpha = 0.5) +
  geom_line(stat = 'smooth',
            method = 'lm',
            linewidth = 1.5,
            alpha = 0.8) +
  ylim(c(0, 3)) +
  theme_bw() +
  ggtitle('Grassland birds')


# negative examples -------------------------------------------------------

tt4 <- dplyr::filter(bird_df2, 
                     avonet.Habitat == 'Forest',
                     !is.na(bird.et.al.GenLength),
                     !is.na(sd_resid_temp))

ggplot(tt4, aes(sd_resid_temp, log(bird.et.al.GenLength),
                col = factor(avonet.Trophic.Level))) +
  geom_point(alpha = 0.1) +
  geom_line(stat = 'smooth',
            method = 'lm',
            linewidth = 1.5,
            alpha = 0.8) +
  ylim(c(0, 3)) +
  theme_bw() +
  ggtitle('Forest birds')


tt5 <- dplyr::filter(bird_df2, 
                     avonet.Habitat == 'Shrubland',
                     !is.na(bird.et.al.GenLength),
                     !is.na(sd_resid_temp))

ggplot(tt5, aes(sd_resid_temp, log(bird.et.al.GenLength),
                col = factor(avonet.Trophic.Level))) +
  geom_point(alpha = 0.3) +
  geom_line(stat = 'smooth',
            method = 'lm',
            linewidth = 1.5,
            alpha = 0.8) +
  ylim(c(0, 3)) +
  theme_bw() +
  ggtitle('Shrubland birds')


# PCA env var -------------------------------------------------------------

bird_df3 <- dplyr::filter(bird_df2, !is.na(mean_temp), 
              !is.na(sd_resid_temp), 
              !is.na(rho_l1_temp),
              avonet.Habitat == 'Forest',
              avonet.Trophic.Level == 'Herbivore')

pr_fit <- prcomp(cbind(bird_df3$mean_temp, 
                       bird_df3$sd_resid_temp,
                       bird_df3$rho_l1_temp), scale = TRUE)

#INTERPRETATION OF PC DEPENDS ON GROUPING
#for forest herbivore
#pos PC1 = low temp, high sd, slightly low rho
#pos PC2 = slightly low sd, low rho
factoextra::fviz_pca_var(pr_fit,
                         #geom = 'arrow',
                         #col.var = "contrib", # Color by contributions to the PC
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = FALSE)
#first 2 PC
bird_df3$PC1 <- pr_fit$x[,1]
bird_df3$PC2 <- pr_fit$x[,2]

summary(lm(log(bird.et.al.GenLength) ~ PC1 + PC2, data = bird_df3))
summary(lm(log(bird.et.al.GenLength) ~ sd_resid_temp, data = bird_df3))

