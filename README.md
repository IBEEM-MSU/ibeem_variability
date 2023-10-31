# IBEEM_variability

## Project Links:

[Project Google Drive](https://drive.google.com/drive/u/1/folders/1hnJP1CRZQSph2M2cIhCEwvujxFOhfM0k)

[Meeting Notes](https://docs.google.com/document/d/1ix5mSCbO7rRCj2juQN04UeabejayMgcYT0-yo08uSZ0/edit)

[Data Documentation](https://docs.google.com/document/d/13Hn0klwabOQzdfCg1W-FIF00j1t_YgT40v-99Fhkt3Y/edit)

[Proposal](https://drive.google.com/file/d/1K0jelDSM3ZlyHDlNI3SRYAbzlrS0hpiU/view?usp=share_link)

[Annotated Bib](https://docs.google.com/document/d/1fn3lzU2IUtvY7Go_O7An_6JLgPk7Zc-OKLlNQYgOq48/edit)


## Repo structure:
* `Scripts/`
  * `1-clean-data/`
    * `1a-clean-Bird-et-al.R` - clean and combine gen time data from Bird et al. 2020
    * `1a-clean-terrestrial-mammal.R` - clean mammal range data 
    * `1b-process-ncdf.R` - create time series of env data, yearly averages for specified months, and sd across specified months (seasonality)
    * `1b-process-ncdf.slurm` - slurm script to submit job
    * `1c-process-ncdf-monthly.R` - create monthly time series of env data, then calculate spectral exponent
    * `1c-create-batch-scripts.sh` - script to create slurm scripts
    * `1c-master-submit.sh` - script to submit batch of slurm scripts
    * `chunks/` - slurm scripts created by `1c-create-batch-scripts.sh`
    * `1d-taxonomic-harmonization.R` - sort out naming differences and save individual range maps for birds
    * `1d-taxonomic-harmonization-terrestrial-mammal.R` - sort out naming differences and save individual range maps for mammals
  * `2-env-metrics/`
    * `2a-env-metrics.R` - calculate env variability metrics
    * `2b-env-metrics-GAM.R` - calculate env variability metrics using GAM detrend
    * `2c-DHI.R` - process DHI data
    * `2d-env-merge.R` - merge env metrics
  * `3-extract-species-env/` - extract env var data from species range
    * `3a-split-sp-ids.R` - generate sets of bird ids for parallel processing
    * `3a-split-sp-ids-terrestrial-mammal.R` - generate sets of mammal ids for parallel processing
    * `3b-extract-avg-within-range.R` - extract environmental covariates from species ranges
    * `3b-extract-avg-within-range-terrestrial-mammal.R` - extract environmental covariates from species ranges
    * `3b-extract-compile.sh` - bash script to iterate through all pieces of data and run `3b-extract.sbatch`
    * `3b-extract.sbatch` - bash script to load one piece of data and run `3b-extract-avg-within-range.R` on HPCC
    * `3b-extract-compile-mammal.sh` - bash script to iterate through all pieces of data and run `3b-extract-mammal.sbatch`
    * `3b-extract-mammal.sbatch` - bash script to load one piece of data and run `3b-extract-avg-within-range-terrestrial-mammal.R` on HPCC
    * `3c-get-master-file-terrestrial-mammal.R` - generate master file with rows for species and columns for environmental and life history data
    * `3c-get-master-file-birdtree.R` - generate master file with rows for species and columns for environmental and life history data. Also pulls in IUCN data.
    + `3c-get-master-file-birdtree.sbatch` - bash script to run `3c-get-master-file-birdtree.R`, since pulling the IUCN data takes a few hours.
  * `4-phylo-cor` - get phylogenetic correlation matrix
    * `4a-bird-get-con-tree-pieces.R` - calculate consensus trees for 10 chunks of 100 phylogenies for birds
    * `4b-bird-get-cor-matrix.R` - calculate final consensus tree and phylo cor matrix for birds 
    * `4a-mammal-get-con-tree-pieces.R` - calculate consensus trees for 10 chunks of 100 phylogenies for mammals
    * `4b-mammal-get-cor-matrix.R` - calculate final consensus tree and phylo cor matrix for mammals 
  * `5-model/` - * = working scripts
    * `Archive/` - old model scripts
      * `5-clutch-oe.R` - imputed clutch size ~ var with observation error
      * `5-gl-pc.R` - gen length ~ PC var
      * `5-gl-varfam.R` - gen length ~ var varying intercepts, varying slopes by family  
      * `5-surv-phylo.R` - imputed survival ~ var with phylo
      * `5-surv-pc-oe.R` - imputed survival ~ PC var with observation error
      * `5-surv-oe.R` - imputed survival ~ var with observation error
      * `5-surv-pc-phylo.R` -imputed surv ~ PC var with phylo
      * `5-surv-phylo-oe.R` - imputed survival ~ var with phylo and observation error  
      * `5-surv-phylo-oe-rs.R` - as above but check benchmarking reduce_sum
      * `5-surv-phylo-oe.slurm` - to run `5-surv-phylo-oe.R`
      * `5-gl-phylo.R` - gen length ~ var + phylo
      * `5-gl-varniche.R` - gen length ~ var + varying intercepts, varying slopes by niche
      * `5-gl.R` - gen length ~ var
      * `5-surv-phylo-vint.R` - survival ~ var + varying int + phylo kappa
      * `5-gl-intvar.R` - gen length ~ var - varying int (by fam/niche) model - not that useful since it's ignoring var across fams/niches
      * `5-gl-mam.R` - gen length ~ var for mammals
      * `bird_df3.rds` - data to avoid running first part of `5-...R` scripts
      * `brms.R` - benchmarking brms models with rstan models
      * `5-gl-phylo-vint.R` - gen length ~ var varying intercepts, phylo kappa
      * `5-ml-phylo-vint.R` - max long ~ var + varying int + phylo kappa
      * `5-ab-phylo-vint.R` - age first breeding ~ var varying intercepts, phylo kappa
      * `5-rl-phylo-vint.R` - rel repro lifespan ~ var + varying int + phylo kappa
      * `5-gl-phylo-vint-gamma.R` - gen length ~ var varying intercepts, phylo kappa with gamma instead of normal
    * `5-ab-phylo-vint-oe.R`* - age first breeding ~ var varying intercepts, phylo kappa; modeled (from Bird et al.) rather than imputed values as in `Archive`
    * `5-ml-phylo-vint-oe.R`* - max long ~ var varying intercepts, phylo kappa; modeled (from Bird et al.) rather than imputed values as in `Archive`
    * `5-s-phylo-vint-oe.R`* - survival ~ var varying intercepts, phylo kappa; modeled (from Bird et al.) rather than imputed values as in `Archive`
    * `5-gl-phylo-vint-mam.R`* - gen length ~ var varying intercepts, phylo kappa mammals
    * `5-iucn-threat.R`* - IUCN threat status ~ variation
  * `6-rasterize` - rasterize ranges
    * `6a-create-rasters.R` - create .tif files birds
    * `6a-raster-comile.sh` - bash script to iterate through all pieces of data and run `6a-raster.sbatch`
    * `6a-raster.sbatch` - bash script to load one piece of data and run `6a-create-rasters.R` on HPCC
    * `6b-stack-rasters.R` - script to stack rasters and produce main .tif
  * `7-figures` - figures
  * `X-explore/`
    * `X-explore-bird.R` - exploratory analyses birds
    * `X-explore-mammal.R` - exploratory analyses mammals
    * `X-explore-env-var.R` - explore environmental variability
  * `Model_files` - Stan model files
* `Archive/` - unused scripts  
* `Data/` (ignored)
  * `L0/` - raw data
    * `DHI/` - Dynamic Habitat Index data
    * `climate/era5/` - raw ERA5 reanalysis data
    * `ranges/BOTW.gdb` - BirdLife range maps
    * `trait/` - raw trait data
  * `L1/`
    * `DHI/` - processed DHI data (cv year and cv season)
    * `climate/era5/` - ERA data averaged over specified months (one value per cell/year)
      * `ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv` - yearly average over all months and sd across months (seasonality) for temp and (sqrt root transform of) precip
      * `Env-spectral-exp-monthly.csv` - spectral exponent of monthly env variables
    * `range/` - bird ranges
    * `range-mammal/` - mammal ranges
    * `trait/` - processed bird traits
    * `trait-mammal/` - processed mammal traits
  * `L2/`
    * `climate/era5/` - env variability metrics per cell (including seasonality)
      * `Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv` - variability for all months
      * `Env-main.csv` - THIS IS THE MAIN MERGED ENV DATA TO USE DOWNSTREAM
    * `range-env-pieces/` - ??? birds
    * `range-env-pieces-mammal/` - ??? mammals
    * `range-raster/` - rasterized ranges (tifs) with gen length and delta haldane (relative temp change) as layers
  * `L3/` - main data files
    * `main-bird-data.csv` - merged bird data
    * `main-mammal-data.csv` - merged mammal data
    * `raster-gl-dh-nsp.tif` - raster of median and sd of gen length and delta haldane (relative temp change), as well as number of species in each cell (5 layers)
    * `bird-phylo-cor-matrix.rda` - bird phylogenetic correlation matrix
    + `bird-consensus-tree.rda` - bird consensus tree
    * `mammal-phylo-cor-matrix.rda` - mammal phylogenetic correlation matrix
    + `mammal-consensus-tree.rda` - mammal consensus tree


## Data sources:
* [ERA5 reanalysis data (temp and precipitation)](https://rda.ucar.edu/datasets/ds633.1/)
* [Bird Life International range maps](http://datazone.birdlife.org/species/requestdis)
* [IUCN mammal range amps](https://www.iucnredlist.org/resources/spatial-data-download)
* [Bird et al. 2020 *Conservation Biology* bird generation time](https://conbio.onlinelibrary.wiley.com/doi/10.1111/cobi.13486)
* [Pacifici et al. 2013 *Nature Conservation* mammal generation time](https://natureconservation.pensoft.net/article/1343/download/pdf/)
* [Elton Traits (Wilman et al. 2014) for birds and mammals](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1917.1)
* [PanTHERIA (Jones et al. 2009) for mammals](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/08-1494.1)
* [AVONET (Tobias et al. 2022) for birds](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898)
* [Cooney et al. 2020 for birds](https://www.nature.com/articles/s41467-020-16257-x)

## Env time series (avg over specified months and sd across months [seasonality]):
* `sbatch variability/Scripts/1-clean-data/1b-process-ncdf.slurm`

## Spectral exponent monthly time series:
* `./variability/Scripts/1-clean-data/1c-master-submit.sh`

## Env variability metics:
* `sbatch variability/Scripts/2-env-metrics/2a-env-1_12.slurm`

## Extracting env and trait data for each species
+ Code for initial processing for birds (breeding season) is in `Scripts/3-extract-species-env`
+ Currently using species with ITIS ids, did not resolve NAs. 
+ Generates a main file with each row a species (with associated ITIS ID and accepted name). Each column contains either the average of one of the environmental variables across the range of the species, or some trait value. 
+ Main file is stored in `/mnt/research/ibeem/data/L2/main-bird-data.csv` and is also in the [`L2` folder on Google drive](https://drive.google.com/drive/folders/1c4dF8AEgOf7zvVUd5GTUT7vzj5hoLuGj). 
+ Trait data in the master file currently include values from Bird et al and AVONET. Can add in ELTON later on if there is something specific we want from there. The column name has the data source as a prefix. See metadata for the data sources for specific info.

