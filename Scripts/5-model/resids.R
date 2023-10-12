range_size_km2

str(bird_df)
bdf <- dplyr::select(bird_df,
              Birdtree_name,
              lGL,
              Habitat, Habitat.Density, 
              Trophic.Level, Trophic.Niche, 
              Tail.Length, Primary.Lifestyle, 
              Mass, Wing.Length, range_size_km2) %>%
  dplyr::mutate(species = stringr::str_to_sentence(gsub(' ', '_', Birdtree_name)))

str(bird_df5)
str(bdf)

t2 <- dplyr::left_join(bird_df5, bdf, by = 'species')
obs_data <- t2[DATA$obs_idx,]
imp_data <- t2[DATA$imp_idx,]

obs_data$resid <- resid_obs
imp_data$resid <- resid_imp

dd <- rbind(obs_data, imp_data)

plot(factor(obs_data$Habitat), obs_data$resid)
plot(factor(obs_data$Habitat.Density), obs_data$resid)
plot(factor(obs_data$Trophic.Level), obs_data$resid)
plot(factor(obs_data$Trophic.Niche), obs_data$resid) #*
plot(obs_data$Tail.Length, obs_data$resid)
plot(factor(obs_data$Primary.Lifestyle), obs_data$resid)
plot(log(obs_data$Wing.Length), obs_data$resid)
plot(log(obs_data$range_size_km2), obs_data$resid)

ggplot(dd, aes(Trophic.Niche, resid)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5)) +
  xlab('') +
  ylab('Residual') +
  ggtitle('Survival ~ Env + Phylo')
  

jpeg('~/Desktop/niche.jpeg', width = 500, height = 500)
plot(factor(dd$Trophic.Niche), dd$resid, las = 2) #*
dev.off()

str(dd)
dplyr::group_by(dd, Trophic.Niche) %>%
  dplyr::summarize(mr = mean(resid),
                   m_t_season = mean(temp_sd_season),
                   m_t_year = mean(temp_sd_year),
                   m_p_season = mean(precip_cv_season),
                   m_p_year = mean(precip_cv_year))
str(dd)

ggplot(dd, aes(temp_sd_season, Phylo_survival, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

ggplot(dd, aes(temp_sd_year, Phylo_survival, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

ggplot(dd, aes(temp_sp_color_month, lGL, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

ggplot(dd, aes(precip_cv_season, Phylo_log_clutch_size, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

ggplot(dd, aes(precip_cv_year, Phylo_log_clutch_size, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

ggplot(dd, aes(precip_sp_color_month, lGL, color = Trophic.Niche)) +
  geom_point(alpha = 0.08) +
  geom_line(stat = 'smooth', method = 'lm',
            linewidth = 2, alpha = 0.5) +
  theme_bw()

unique(dd$Trophic.Niche)
dplyr::filter(dd, is.na(Trophic.Niche))$species

      
plot(factor(obs_data$Trophic.Niche), obs_data$resid, ylim = c(-3, 3)) #*
plot(factor(dd$Trophic.Niche), dd$resid, ylim = c(-3, 3)) #*
