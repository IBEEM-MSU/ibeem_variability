# TITLE:            Calculate uncertainty in gen length
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       survival, breeding, longevity data (Bird et al. 2020)
# DATA OUTPUT:      R^2 for gen length
# DATE:             July 2024
# OVERVIEW:         Calculate uncertainty for generation length values given by Bird et al. 2020 using simulation


# load environment variables ------------------------------------------------

source('./Scripts/0-config.R')


# load packages -----------------------------------------------------------

library(tidyverse)


# load bird data ---------------------------------------------------------------

#Bird et al. Supp Table 4
#Bird, J. P., R. Martin, H. R. Akçakaya, J. Gilroy, I. J. Burfield, S. Garnett, A. Symes, J. Taylor, Ç. H. Şekercioğlu, and S. H. M. Butchart. 2020. Generation lengths of the world’s birds and their implications for extinction risk. Conservation Biology. 34(5), 1252-1261. (DOI: doi.org/10.1111/cobi.13486)

table4 <- readxl::read_excel(paste0(dir, 'data/L0/Bird_et_al_2020/cobi13486-sup-0004-tables4.xlsx'))

t4_mod <- dplyr::rename(table4, 
                        Final_survival = 'Adult survival',
                        Final_age_first_breeding = 'Age at first breeding',
                        Final_max_longevity = 'Maximum longevity',
                        Sci_name = 'Scientific name') %>%
  dplyr::select(Order, Family, Genus, Sci_name, Final_survival, 
                Final_age_first_breeding, Final_max_longevity, GenLength)


# function from Bird et al. Eq. 1 ----------------------------------------------

#calc gen length from Fi, L, and S
g_fun <- function(Fi, L, S)
{
  #numerator
  num <- 0
  #denomonator
  den <- 0
  for (i in Fi:L)
  {
    #survival to year i
    lx <- S^(i - 1)
    num <- num + (lx * i)
    den <- den + lx
  }
  g <- num/den
  return(g)
}


# simed or obs data -------------------------------------------------------

afb <- t4_mod$Final_age_first_breeding
ml <- t4_mod$Final_max_longevity
surv <- t4_mod$Final_survival


# R2 from Bird et al. Table 2 ---------------------------------------------

# #r2 = 1 - (resid var / total var)
# #resid var / total var = (1 - r2)
# #(1 - r2) * total var = resid var
sd_Fi <- sqrt((1 - 0.89) * var(afb))
sd_L <- sqrt((1 - 0.66) * var(ml))
sd_S <- sqrt((1 - 0.77) * var(surv))


# 'truth' -----------------------------------------------------------------

tg <- rep(NA, length(surv))
for (j in 1:length(surv))
{
  #j <- 1
  tg[j] <- g_fun(Fi = afb[j], 
                 L = ml[j], 
                 S = surv[j])
}


# uncertainty -------------------------------------------------------------

N <- 100
#rows are 'species', col are replicates
gmat <- matrix(NA, nrow = length(surv), ncol = N)
r2 <- rep(NA, N)
for (i in 1:N)
{
  #i <- 1
  print(paste0('processing ', i, ' of ', N))
  #for each 'species'
  for (j in 1:length(surv))
  {
    #j <- 1
    #random draw for uncertainty
    Fi_draw <- rnorm(1, afb[j], sd_Fi)
    L_draw <- rnorm(1, ml[j], sd_L)
    #beta for S, normal for Fi and L
    # https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
    b_alpha <- (((1 - surv[j]) / sd_S^2) - (1 / surv[j])) * surv[j]^2
    b_beta <- b_alpha * ((1 / surv[j]) - 1)
    S_draw <- rbeta(1, b_alpha, b_beta)
    #get gl
    gmat[j, i] <- g_fun(Fi = Fi_draw, L = L_draw, S = S_draw)
  }
  #calc r^2 for this rep
  ff <- summary(lm(gmat[,i] ~ tg))
  r2[i] <- ff$r.squared
}


# r2 --------------------------------------------------------------

#0.87
median(r2)

