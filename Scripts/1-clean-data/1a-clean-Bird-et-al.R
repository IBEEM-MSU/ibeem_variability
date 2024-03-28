# TITLE:            Clean Bird et al. data 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Kelly Kapsar, Phoebe L. Zarnetske
# DATA INPUT:       Excel tables with bird gen time from Bird et al. 2020 Conservation Biology
# DATA OUTPUT:      Clean csv with bird generation length data 
# DATE:             May 2023 
# OVERVIEW:         Combine/clean bird gen time data 


# specify dir -------------------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data ------------------------------------------------------------

#from Bird et al. 2020: 
#Published values of F, L and S and the original source (Appendix S1), published (observed), modelled and extrapolated values of F, L and S (Appendix S2), a database of avian life-history traits recorded for individual species (Appendix S3), calculated GFLS, GFL, GFS and G for all birds (noting those considered priorities for research prior to the adoption of the G for Red List assessments) (Appendix S4), and revised Red List trend estimates with predicted qualifying Red List categories for all birds (Appendix S5) are available online. The authors are solely responsible for the content and functionality of these materials. Queries (other than absence of the material) should be directed to the corresponding author.

#ignored Data dir on git...small size though so could add it if desired
#recorded
table3 <- readxl::read_excel(paste0(dir, 'Data/L0/Bird_et_al_2020/cobi13486-sup-0003-tables3.xlsx'))
#modeled - should use GenLength
table4 <- readxl::read_excel(paste0(dir, 'Data/L0/Bird_et_al_2020/cobi13486-sup-0004-tables4.xlsx'))


# combine -----------------------------------------------------------------

#notes: some species name in table3 are not correct... but species order appears to be the same
NROW(table3) == NROW(table4)

#rename cols and keep only necessary ones
t3_mod <- dplyr::rename(table3, 
                        Measured_clutch_size = 'Mean clutch size',
                        Measured_survival = Survival,
                        Measured_age_first_breeding = 'Age at first breeding',
                        Measured_max_longevity = 'Max longevity') %>%
  dplyr::select(Measured_clutch_size, Measured_survival, 
                Measured_age_first_breeding, Measured_max_longevity)

t4_mod <- dplyr::rename(table4, 
                        Modeled_survival = 'Adult survival',
                        Modeled_age_first_breeding = 'Age at first breeding',
                        Modeled_max_longevity = 'Maximum longevity',
                        Sci_name = 'Scientific name') %>%
  dplyr::select(Order, Family, Genus, Sci_name, Modeled_survival, 
                Modeled_age_first_breeding, Modeled_max_longevity, GenLength)

#combine
comb <- cbind(t4_mod, t3_mod) %>%
  dplyr::relocate(c(Measured_survival, 
                    Measured_age_first_breeding,
                    Measured_max_longevity), .before = Modeled_survival)


# write out ---------------------------------------------------------------

write.csv(comb, paste0(dir, 'Data/L0/trait/Bird_et_al_gen_length_birds.csv'))

          