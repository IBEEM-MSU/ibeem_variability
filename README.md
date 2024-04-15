# Variability_life_history

Code for Youngflesh et al. *In prep*

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species.

**Associated publications:**

TBD


**Project Links:**

[Project Google Drive](https://drive.google.com/drive/u/1/folders/1hnJP1CRZQSph2M2cIhCEwvujxFOhfM0k)

[Meeting Notes](https://docs.google.com/document/d/1ix5mSCbO7rRCj2juQN04UeabejayMgcYT0-yo08uSZ0/edit)

[Data Documentation](https://docs.google.com/document/d/13Hn0klwabOQzdfCg1W-FIF00j1t_YgT40v-99Fhkt3Y/edit)

[Proposal](https://drive.google.com/file/d/1K0jelDSM3ZlyHDlNI3SRYAbzlrS0hpiU/view?usp=share_link)

[Annotated Bib](https://docs.google.com/document/d/1fn3lzU2IUtvY7Go_O7An_6JLgPk7Zc-OKLlNQYgOq48/edit)


**Repository structure:**
* `Scripts/`
  * `1-clean-data/`
    * `1a-clean-Bird-et-al.R` - clean and combine gen time data from Bird et al. 2020
    * `1b-process-ncdf.R` - create time series of env data, yearly averages for specified months, and sd across specified months (seasonality)
    * `1b-process-ncdf.slurm` - slurm script to submit job
    * `1c-taxonomic-harmonization.R` - sort out naming differences and save individual range maps for birds
  * `2-env-metrics/`
    * `2a-env-metrics.R` - calculate env variability metrics
    * `2b-env-merge.R` - merge env metrics
  * `3-extract-species-env/` - extract env var data from species range
    * `3a-split-sp-ids.R` - generate sets of bird ids for parallel processing
    * `3b-extract-avg-within-range.R` - extract environmental covariates from species ranges
    * `3b-extract-compile.sh` - bash script to iterate through all pieces of data and run `3b-extract.sbatch`
    * `3b-extract.sbatch` - bash script to load one piece of data and run `3b-extract-avg-within-range.R` on HPCC HPCC
    * `3c-get-master-file-birdtree.R` - generate master file with rows for species and columns for environmental and life history data. Also pulls in IUCN data.
    + `3c-get-master-file-birdtree.sbatch` - bash script to run `3c-get-master-file-birdtree.R`, since pulling the IUCN data takes a few hours.
  * `4-phylo-cor` - get phylogenetic correlation matrix
    * `4a-bird-get-con-tree-pieces.R` - calculate consensus trees for 10 chunks of 100 phylogenies for birds
    * `4b-bird-get-cor-matrix.R` - calculate final consensus tree and phylo cor matrix for birds 
  * `5-model/` -
    * `5-gl-phylo-vint.R` - fit gl ~ env var
    * `5-gl-phylo-vint.slurm` - slurm script to submit job
  * `6-rasterize` - rasterize ranges
    * `6a-create-rasters.R` - create .tif files birds
    * `6a-raster-comile.sh` - bash script to iterate through all pieces of data and run `6a-raster.sbatch`
    * `6a-raster.sbatch` - bash script to load one piece of data and run `6a-create-rasters.R` on HPCC
    * `6b-stack-rasters.R` - script to stack rasters and produce main .tif
  * `7-figures` - figures
  * `Model_files` - Stan model files
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
    * `trait/` - processed bird traits
  * `L2/`
    * `climate/era5/` - env variability metrics per cell (including seasonality)
      * `Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv` - variability for all months
      * `Env-main.csv` - THIS IS THE MAIN MERGED ENV DATA TO USE DOWNSTREAM
    * `range-env-pieces/` - ??? birds
    * `range-env-pieces-mammal/` - ??? mammals
    * `range-raster/` - rasterized ranges (tifs) with gen length and delta haldane (relative temp change) as layers
  * `L3/` - main data files
    * `main-bird-data.csv` - merged bird data
    * `raster-gl-dT-dP-nsp.tif` - raster of median and sd of gen length and delta haldane (relative temp change), as well as number of species in each cell (5 layers)
    * `bird-phylo-cor-matrix.rda` - bird phylogenetic correlation matrix
    * `bird-consensus-tree.rda` - bird consensus tree
