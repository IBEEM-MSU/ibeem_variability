# TITLE:            Figure S2: Pairwise niche figure   
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Kelly Kapsar, Phoebe L. Zarnetske
# DATA INPUT:        
# DATA OUTPUT:      
# DATE:             January 2024 
# OVERVIEW:         

# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# dir <- '/mnt/research/ibeem/variability/'
# sc_dir <- '/mnt/home/ccy/variability/'
gl_run_date <- '2023-10-17'


# read in data ------------------------------------------------------------

#pairwise difference from model fit
pw_df2 <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
                          '/pairwise_gamma.rds'))


#plot of mean percent difference
p1 <- ggplot(pw_df2, aes(niche_names.x, niche_names.y, fill = abs(mean))) +
  geom_tile(color = 'black', linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(vjust = 1, size = 12, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())

ggsave(p1, filename = '~/Desktop/niche_gen_length.pdf')

#plot of prob greater than 0 for estimate - SHOULD BE DIVERGING COLOR SCALE
p2 <- ggplot(pw_df2, aes(niche_names.x, niche_names.y, fill = Pg0)) +
  geom_tile(color = 'black', linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(vjust = 1, size = 12, hjust = 1)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())
