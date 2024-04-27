# Variability_life_history

**DOI HERE**

Code for Youngflesh et al. *In review*

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species.

**Associated publications:**

Youngflesh, C, K Kapsar, A Uscanga, PJ Williams, JW Doser, L Kounta, PL Zarnetske. Environmental variation shapes life history of the world’s birds. *In Review*


**Repository structure:**
* `Scripts/`
  * `0-config.R` - congfiguration file specifying location of data directory and model run date
  * `0-csv-to-tif.R` - convert csv files back into geospatial format (tif) for analysis  
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
    * `3b-extract.sbatch` - bash script to load one piece of data and run `3b-extract-avg-within-range.R` on HPCC
    * `3c-get-master-file-birdtree.R` - generate master file with rows for species and columns for environmental and life history data
    + `3c-get-master-file-birdtree.sbatch` - bash script to run `3c-get-master-file-birdtree.R`
  * `4-phylo-cor` - get phylogenetic correlation matrix
    * `4a-bird-get-con-tree-pieces.R` - calculate consensus trees for 10 chunks of 100 phylogenies for birds
    * `4b-bird-get-cor-matrix.R` - calculate final consensus tree and phylo cor matrix for birds 
  * `5-model/`
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
    * `climate/era5/` - raw ERA5 reanalysis data
    * `ranges/BOTW.gdb` - BirdLife range maps
    * `trait/` - raw trait data
  * `L1/`
    * `climate/era5/` - ERA data averaged over specified months (one value per cell/year)
      * `ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv` - yearly average over all months and sd across months (seasonality) for temp and (sqrt root transform of) precip
    * `range/` - bird ranges
    * `trait/` - processed bird traits
  * `L2/`
    * `climate/era5/` - env variability metrics per cell
      * `Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv` - variability for all months
      * `Env-main.csv` - main env data
    * `range-env-pieces/` - IDs for species
    * `range-raster/` - rasterized ranges (tifs) with gen length and delta metrics (relative abiotic change) as layers
  * `L3/` - main data files
    * `main-bird-data.csv` - merged bird data
    * `raster-gl-dT-dP-nsp.tif` - raster of median and sd of gen length and delta metrics (relative abiotic change), as well as number of species in each cell (5 layers)
    * `bird-phylo-cor-matrix.rda` - bird phylogenetic correlation matrix
    * `bird-consensus-tree.rda` - bird consensus tree
    * `final-bird-data-for-archival.csv` - archival version of `main-bird-data.csv`
    * `env_var.csv` - archival environmental variability grid data
    * `delta.csv` - archival delta_T and delta_P grid data 
