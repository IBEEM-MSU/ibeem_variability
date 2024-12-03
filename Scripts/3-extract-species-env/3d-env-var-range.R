# TITLE:            Environmental variation over species' ranges
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Environmental data
# DATA OUTPUT:      Processed bird/enviromental data
# DATE:             July 2024
# OVERVIEW:         Plot degree of variation across range in env variability


rm(list = ls())


# load environment variables ------------------------------------------------

source('./Scripts/0-config.R')


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data ------------------------------------------------------------

or_excl <- c('Sphenisciformes', #penguins
             'Procellariiformes', #tubenoses
             'Pelecaniformes', #pelicans
             'Suliformes', #gannets/boobies
             'Phaethontiformes', #tropicbirds
             'Charadriiformes')#, #skuas, gulls, terns, skimmers, auks

'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1,
                !is.na(temp_sd_year))


# plot variation for each species -----------------------------------------

df2 <- data.frame(temp_sd_year = bird_df$temp_sd_year,
                  temp_sd_year_LCI = bird_df$temp_sd_year - 
                    bird_df$temp_sd_year_sd_space,
                  temp_sd_year_UCI = bird_df$temp_sd_year + 
                    bird_df$temp_sd_year_sd_space,
                  temp_sd_season = bird_df$temp_sd_season,
                  temp_sd_season_LCI = bird_df$temp_sd_season -
                    bird_df$temp_sd_season_sd_space,
                  temp_sd_season_UCI = bird_df$temp_sd_season +
                    bird_df$temp_sd_season_sd_space,
                  precip_cv_year = bird_df$precip_cv_year,
                  precip_cv_year_LCI = bird_df$precip_cv_year -
                    bird_df$precip_cv_year_sd_space,
                  precip_cv_year_UCI = bird_df$precip_cv_year +
                    bird_df$precip_cv_year_sd_space,
                  precip_cv_season = bird_df$precip_cv_season,
                  precip_cv_season_LCI = bird_df$precip_cv_season -
                    bird_df$precip_cv_season_sd_space,
                  precip_cv_season_UCI = bird_df$precip_cv_season +
                    bird_df$precip_cv_season_sd_space) %>%
  dplyr::mutate(rng_temp_sd_year = temp_sd_year_UCI - 
                  temp_sd_year_LCI,
                rng_temp_sd_season = temp_sd_season_UCI - 
                  temp_sd_season_LCI,
                rng_precip_cv_year = precip_cv_year_UCI - 
                  precip_cv_year_LCI,
                rng_precip_cv_season = precip_cv_season_UCI - 
                  precip_cv_season_LCI)

#temp sd year
t_sd_year_plt <- df2 %>%
  # dplyr::arrange(temp_sd_year) %>%
  dplyr::arrange(rng_temp_sd_year) %>%
  ggplot() +
  geom_point(aes(x = 1:NROW(df2), 
                 temp_sd_year),
             alpha = 0.1) +
  geom_errorbar(aes(x = 1:NROW(df2), 
                    ymin = temp_sd_year_LCI,
                    ymax = temp_sd_year_UCI),
                alpha = 0.05) +
  # ggtitle('temp sd year') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks.y = element_line(linewidth = 1.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank())
# axis.title.y = element_text(size = 18))


#temp sd season
t_sd_season_plt <- df2 %>%
  # dplyr::arrange(temp_sd_season) %>%
  dplyr::arrange(rng_temp_sd_season) %>%
  ggplot() +
  geom_point(aes(x = 1:NROW(df2), 
                 temp_sd_season),
             alpha = 0.1) +
  geom_errorbar(aes(x = 1:NROW(df2), 
                    ymin = temp_sd_season_LCI,
                    ymax = temp_sd_season_UCI),
                alpha = 0.05) +
  # ggtitle('temp sd season') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks.y = element_line(linewidth = 1.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank())
# axis.title.y = element_text(size = 18))

#precip cv year
p_cv_year_plt <- df2 %>%
  # dplyr::arrange(precip_cv_year) %>%
  dplyr::arrange(rng_precip_cv_year) %>%
  ggplot() +
  geom_point(aes(x = 1:NROW(df2), 
                 precip_cv_year),
             alpha = 0.1) +
  geom_errorbar(aes(x = 1:NROW(df2), 
                    ymin = precip_cv_year_LCI,
                    ymax = precip_cv_year_UCI),
                alpha = 0.05) +
  # ggtitle('precip cv year') +
  theme_bw() +
  coord_cartesian(ylim =  c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks.y = element_line(linewidth = 1.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank())
# axis.title.y = element_text(size = 18))

#precip cv season
p_cv_season_plt <- df2 %>%
  # dplyr::arrange(precip_cv_season) %>%
  dplyr::arrange(rng_precip_cv_season) %>%
  ggplot() +
  geom_point(aes(x = 1:NROW(df2), 
                 precip_cv_season),
             alpha = 0.1) +
  geom_errorbar(aes(x = 1:NROW(df2), 
                    ymin = precip_cv_season_LCI,
                    ymax = precip_cv_season_UCI),
                alpha = 0.05) +
  # ggtitle('precip cv season') +
  coord_cartesian(ylim =  c(0,1.75)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks.y = element_line(linewidth = 1.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank())
        # axis.title.y = element_text(size = 18))

ggsave(filename = paste0(dir, '/Results/env-var-rng/t_sd_year.png'),
       t_sd_year_plt,
       height = 5,
       width = 5)
ggsave(filename = paste0(dir, '/Results/env-var-rng/t_sd_season.png'),
       t_sd_season_plt,
       height = 5,
       width = 5)
ggsave(filename = paste0(dir, '/Results/env-var-rng/p_cv_year.png'),
       p_cv_year_plt,
       height = 5,
       width = 5)
ggsave(filename = paste0(dir, '/Results/env-var-rng/p_cv_season.png'),
       p_cv_season_plt,
       height = 5,
       width = 5)



