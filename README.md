# IBEEM_variability

## Project Links:

[Project Google Drive](https://drive.google.com/drive/u/1/folders/1hnJP1CRZQSph2M2cIhCEwvujxFOhfM0k)

[Meeting Notes](https://docs.google.com/document/d/1ix5mSCbO7rRCj2juQN04UeabejayMgcYT0-yo08uSZ0/edit)

[Data Documentation](https://docs.google.com/document/d/13Hn0klwabOQzdfCg1W-FIF00j1t_YgT40v-99Fhkt3Y/edit)

[Proposal](https://drive.google.com/file/d/1K0jelDSM3ZlyHDlNI3SRYAbzlrS0hpiU/view?usp=share_link)


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
    * `3c-get-master-file.R` - generate master file with rows for species and columns for environmental and life history data
  * `4-explore-bird.R` - explore joined bird life history/env data
  * `X-explore-env-var.R` - explore environmental variability metrics
* `Archive/` - ununsed scripts  
* `Data/` (ignored)
  * `L0/` - raw data
    * `DHI/` - Dynamic Habitati Index data
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

