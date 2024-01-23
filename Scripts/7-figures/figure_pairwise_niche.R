####################
# Figure SX - pairwise niche figs
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




#plot of mean percent difference
ggplot(pw_df2, aes(niche_names.x, niche_names.y, fill = abs(mean))) +
  geom_tile(color = 'black', linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(vjust = 1, size = 12, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())

#plot of prob greater than 0 for estimate
ggplot(pw_df2, aes(niche_names.x, niche_names.y, fill = Pg0)) +
  geom_tile(color = 'black', linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(vjust = 1, size = 12, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())