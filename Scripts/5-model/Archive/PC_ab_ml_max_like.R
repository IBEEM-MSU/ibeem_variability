####################
# Max like PGLS - PCs (Ab and Ml) ~ env
####################


# specify dir -------------------------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'
sc_dir <- '/mnt/home/ccy/variability/'
run_date <- '2023-10-24'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)
library(ape)
library(geiger)
library(phytools)


# load bird data ---------------------------------------------------------------

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
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtreeX.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1) %>%
  #filter out precip outlier
  dplyr::filter(precip_cv_season < 2.5) %>%
  dplyr::mutate(lMass = log(Mass),
                lGL = log(GenLength),
                lAb = log(Modeled_age_first_breeding),
                lMl = log(Modeled_max_longevity),
                lCs = log(Measured_clutch_size),
                # lMl = log(Measured_max_longevity),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df, species,
                          lMass,
                          lGL,
                          Modeled_survival,
                          Measured_survival,
                          Measured_max_longevity,
                          Measured_age_first_breeding,
                          lAb,
                          lMl,
                          lCs,
                          Trophic_niche = Trophic.Niche,
                          Order,
                          Family,
                          cen_lon,
                          cen_lat,
                          temp_mean,
                          temp_sd_year,
                          temp_sd_season,
                          precip_cv_year,
                          precip_cv_season)


NROW(dplyr::filter(bird_df, !is.na(Measured_survival))) / NROW(bird_df)
NROW(dplyr::filter(bird_df, !is.na(Measured_max_longevity))) / NROW(bird_df)
NROW(dplyr::filter(bird_df, !is.na(Measured_age_first_breeding))) / NROW(bird_df)

sum(bird_df$Measured_survival > 0.3, na.rm = TRUE) / 
  NROW(dplyr::filter(bird_df, !is.na(Measured_survival)))

hist(bird_df$Measured_survival^bird_df$Measured_max_longevity, 
     breaks = 1000, xlim = c(0, 0.01))

plot(bird_df$Measured_max_longevity, bird_df$Measured_survival^bird_df$Measured_max_longevity,
     xlim = c(0, 40))
hist(bird_df$Measured_max_longevity, breaks = 200)



t3 <- dplyr::mutate(bird_df, 
              lAb = log(Modeled_age_first_breeding),
              lMl = log(Modeled_max_longevity)) %>%
  dplyr::select(Modeled_survival, lAb, lMl)

plot(t3)


# niche levels ------------------------------------------------------------

#assign missing species (owls) to Vertivore
bird_df2$Trophic_niche[which(is.na(bird_df2$Trophic_niche))] <- "Vertivore"

bird_df2$niche_idx <- as.numeric(factor(bird_df2$Trophic_niche))
niche_names <- levels(factor(bird_df2$Trophic_niche))


# clutch size -------------------------------------------------------------

# plot(bird_df2$lCs, bird_df2$lAb)
# plot(bird_df2$lCs, bird_df2$lMl)


# PCA ---------------------------------------------------------------------

# bird_df4 <- dplyr::filter(bird_df2, !is.na(lCs))
bird_df4 <- dplyr::filter(bird_df2, !is.na(Measured_age_first_breeding))
# bird_df4 <- bird_df2

#stronger cor between Ab and Ml than between Lh and Cs
cor(dplyr::select(bird_df4, lMl, lAb, lCs))


#same results as raw, essentially
# tt_pca <- dplyr::select(bird_df4, 
#                         lAb,
#                         lMl) %>%
tt_pca <- dplyr::select(bird_df4,
                          lAb,
                          lMl,
                        lCs) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')


factoextra::fviz_pca_biplot(tt_pca,
                            axes = c(1,2), 
                            label = "var",
                            habillage = bird_df4$Order)#,
                            # addEllipses = TRUE, 
                            #ellipse.level = 0.95)

# +PC1 = longer max long and later first breed
# +PC2 = longer max long and earlier first breed

bird_df4$PC1 <- tt_pca$x[,1]
bird_df4$PC2 <- tt_pca$x[,2]
bird_df4$PC3 <- tt_pca$x[,3]

plot(bird_df4$lAb, bird_df4$lMl,
     xlab = 'log(Age First Breeding)',
     ylab = 'log(Max Longevity)')


# phylo -------------------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#prune tree
#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df4$species)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label), 
                           data.frame(idx = 1:NROW(bird_df4), bird_df4), 
                           by = 'species')

#apply
bird_df5 <- bird_df4[j_idx3$idx,]

#make tree binary
pr_tree2 <- ape::multi2di(pr_tree)

#make response into matrix with species as rownames
dd <- dplyr::select(bird_df5, 
                    lAb) %>%
  as.matrix()
row.names(dd) <- bird_df5$species

#get estimate of Pagel's kappa to scale phylogeny
fit_ka <- geiger::fitContinuous(pr_tree2, dd[,'lAb'], model = "kappa")

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'kappa', 
                               kappa = fit_ka$opt$kappa, sigsq = fit_ka$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)

CM <- nlme::corSymm(Rho[lower.tri(Rho)], fixed = T)

library(nlme)
pgls_fit1 <- nlme::gls(lAb ~ lMass + temp_sd_season + temp_sd_year + precip_cv_season + precip_cv_year,
                       correlation = CM,
                       data = bird_df5,
                       method = "REML")
summary(pgls_fit1)

bird_df5$residA <- residuals(lm(lAb ~ lMass, data = bird_df5))
bird_df5$residM <- residuals(lm(lMl ~ lMass, data = bird_df5))
plot(bird_df5$temp_sd_season, bird_df5$residM)




#from Morrow et al. 2021 Proc B (originally from Charnov pubs)
# LRE (Lifetime Reproductive Effort) represents the energetic trade-off of reproductive effort versus adult mortality
# - clutch size x clutches per year x (mass at independence / adult mass) x adult lifespan
# RRL (Relative Reproductive Lifespan) captures the trade-off in time spent in growth/development versus reproduction
# - adult lifespan / age maturity
# ROS (Relative Offspring Size) is related to trade-offs in the size, number and survivorship of offspring. 
# - mass at independence / adult mass

#RRL 
rrl <- (exp(bird_df5$lMl) - exp(bird_df5$lAb)) / exp(bird_df5$lAb)

rrl <- (bird_df5$Measured_max_longevity - bird_df5$Measured_age_first_breeding) / 
          bird_df5$Measured_age_first_breeding

f1 <- lm(log(rrl) ~ lMass + 
             temp_sd_season + temp_sd_year + 
             precip_cv_season + precip_cv_year,
           data = bird_df5)
summary(f1)
#greater env var, longer repro lifespan
car::crPlots(f1)

#similar to LRE from Morrow et al. 2021 Proc B
lre <- exp(bird_df5$lCs) * (exp(bird_df5$lMl) - exp(bird_df5$lAb))

f2 <- lm(log(lre) ~ lMass + 
           temp_sd_season + temp_sd_year + 
           precip_cv_season + precip_cv_year,
         data = bird_df5)
summary(f2)
plot(bird_df5$lMass, log(lre))
# More var associated with greater max long and later age first breeding

# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)      -2.6309402 0.3514329 -7.48632  0.0000
# lMass             0.7335354 0.0128081 57.27134  0.0000
# temp_sd_season    0.0198306 0.0038474  5.15432  0.0000
# temp_sd_year      0.0222786 0.0586424  0.37991  0.7040
# precip_cv_season -0.0258687 0.0335943 -0.77003  0.4413
# precip_cv_year    0.2927592 0.1492433  1.96162  0.0498

str(bird_df5)
bird_na <- dplyr::filter(bird_df5, Order == 'Passeriformes',
              cen_lon < -60, cen_lon > -130,
              cen_lat > 30)

bird_sa <- dplyr::filter(bird_df5, Order == 'Passeriformes',
                         cen_lon < -60, cen_lon > -130,
                         cen_lat > -20, cen_lat < 20)

f <- lm(bird_na$lGL ~ bird_na$lMass)
naf <- residuals(f)
f2 <- lm(bird_sa$lGL ~ bird_sa$lMass)
saf <- residuals(f2)

ml_f <- lm(bird_na$lMl ~ bird_na$lMass)
ml_naf <- residuals(ml_f)
ml_f2 <- lm(bird_sa$lMl ~ bird_sa$lMass)
ml_saf <- residuals(ml_f2)

ab_f <- lm(bird_na$lAb ~ bird_na$lMass)
ab_naf <- residuals(ab_f)
ab_f2 <- lm(bird_sa$lAb ~ bird_sa$lMass)
ab_saf <- residuals(ab_f2)

plot(density(bird_na$lMass))
lines(density(bird_sa$lMass), col = 'red')

plot(density(naf))
lines(density(saf), col = 'red')

plot(density(ml_naf))
lines(density(ml_saf), col = 'red')

plot(density(ab_naf))
lines(density(ab_saf), col = 'red')

t1 <- lm(lMl ~ lMass + temp_sd_season + temp_sd_year, data = bird_na)
t2 <- lm(lMl ~ lMass + temp_sd_season + temp_sd_year, data = bird_sa)
summary(lm(lCs ~ abs(cen_lat) + lMass, data = bird_df5))
plot(bird_df5$cen_lat, bird_df5$lCs)
plot(bird_na$temp_sd_year, ab_naf)
plot(density(bird_na$temp_sd_season))
lines(density(bird_sa$temp_sd_season), col = 'red')

plot(density(bird_na$temp_sd_year))
lines(density(bird_sa$temp_sd_year), col = 'red')


# PC2

#make response into matrix with species as rownames
dd <- dplyr::select(bird_df5, 
                    PC2) %>%
  as.matrix()
row.names(dd) <- bird_df5$species

#get estimate of Pagel's kappa to scale phylogeny
fit_ka <- geiger::fitContinuous(pr_tree2, dd[,'PC2'], model = "kappa")

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'kappa', 
                               kappa = fit_ka$opt$kappa, sigsq = fit_ka$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)

CM <- nlme::corSymm(Rho[lower.tri(Rho)], fixed = T)

library(nlme)
pgls_fit2 <- nlme::gls(PC2 ~ lMass + temp_sd_season + temp_sd_year + precip_cv_season + precip_cv_year,
                      correlation = CM,
                      data = bird_df5,
                      method = "REML")
summary(pgls_fit2)

# More var associated with greater max long and earlier age first breeding

# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)       0.26601462 0.24218538  1.098393  0.2721
# lMass            -0.03783045 0.00929741 -4.068922  0.0000
# temp_sd_season    0.03974554 0.00282129 14.087706  0.0000
# temp_sd_year      0.17211086 0.04307081  3.995997  0.0001
# precip_cv_season -0.00709479 0.02462611 -0.288100  0.7733
# precip_cv_year   -0.18801015 0.10924910 -1.720931  0.0853


