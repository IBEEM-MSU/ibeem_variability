# Variability_life_history

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species. 

## Table of Contents
- [Overview](#Overview)
- [Funding](#Funding)
- [Associated publications](#Associated-publications)
- [Location of data](#Location-of-data)
- [Spatiotemporal extent and resolution](#Spatiotemporal-extent-and-resolution)
- [Workflow](#Workflow)
- [File naming conventions](#File-naming-conventions)
- [Usage](#Usage)
- [License](#License)
- [Contributors](#Contributors)
- [Contact Information](#Contact-information)

## Overview

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species. The project is part of the [Institute for Biodiversity, Ecology, Evolution, and Macrosystems (IBEEM)](https://ibeem.msu.edu/) at Michigan State University, and led by IBEEM Postdoctoral Fellows. The research uses data from open access biodiversity repositories to evaluate how patterns of environmental variability within a year (i.e., seasons) and between years affect species' life histories. The analysis investigates environmental variability within and between years across the ranges of 7,477 non-migratory and non-marine avian species in order to evaluate the impacts of temporal environmental variability on species' pace of life. Environmental data are from the ERA5 climate reanalysis, bird trait and distribution data are from AVONET, BirdTree, and BirdLife databases as well as previously published research. Environmental variability values were extracted across individual species' ranges in order to quantify how environmental change has affected avian species' pace of life around the world. For more on the project see the [project research summary](https://ibeem.msu.edu/research--news.html) and the publications linked below. 

## Funding

Funding is provided by a Michigan State University Strategic Partnership Grant for the [Institute for Biodiversity, Ecology, Evolution, and Macrosystems (IBEEM)](https://ibeem.msu.edu/).

## Associated publications

Youngflesh, C, K Kapsar, A Uscanga, PJ Williams, JW Doser, L Kounta, PL Zarnetske. Environmental variation shapes life history of the world’s birds. *In Review*

## Location of data

Youngflesh, C., K. Kapsar, A. Uscanga, P.J. Williams, J.W. Doser, L. Kounta, and P.L. Zarnetske. 2024. Inter- and intra-annual temperature and precipitation variability (1950-2022) across the ranges of non-migratory birds and their association with generation length ver 2. Environmental Data Initiative. [https://doi.org/10.6073/pasta/fc057e288924cb09fc3eb5092856d99c](https://doi.org/10.6073/pasta/fc057e288924cb09fc3eb5092856d99c) (Accessed 2024-09-11).

## Spatiotemporal extent and resolution

- Spatial extent: global
- Spatial resolution: 0.25° (climate data)
- Temporal extent: 1950–2022 (climate data)

## Usage

Software used to analyze data includes: R version 4.3.2

## Workflow

The workflow for this repository follows the guidelines set out by the Environmental Data Initiative (EDI). Briefly, this involves aligning with FAIR data practices, and employing a workflow that uses different levels for harmonization and derived data products. Data are read in as raw data at Level 0 (L0). Level 1 (L1) data represent cleaned L0 data. Level 2 data are data merged or otherwise derived from two or more L1 data, etc.

**Repository structure:**

*Note: Scripts are designed to be run in sequence from 0-XXXX.R -> 1-XXXX.R -> ...*

* `Scripts/`
  * `0-config.R` - congfiguration file specifying location of data directory and model run date
  * `0-csv-to-tif.R` - convert csv files back into geospatial format (tif) for analysis  
  * `1-clean-data/`
    * `1a-clean-Bird-et-al.R` - clean and combine gen time data from Bird et al. 2020
    * `1b-process-ncdf.R` - create time series of env data, yearly averages for specified months, and sd across specified months (seasonality)
    * `1b-process-ncdf.slurm` - slurm script to submit job
    * `1c-taxonomic-harmonization.R` - sort out naming differences and save individual range maps for birds
    * `1d-uncertainty-gl.R` - calculate gen length uncertainty
  * `2-env-metrics/`
    * `2a-env-metrics.R` - calculate env variability metrics
    * `2b-env-merge.R` - merge env metrics
    * `2c-env-var-time.R` - explore trend in env variation over time
  * `3-extract-species-env/` - extract env var data from species range
    * `3a-split-sp-ids-birdtree.R` - generate sets of bird ids for parallel processing
    * `3b-extract-avg-within-range-birdtree.R` - extract environmental covariates from species ranges
    * `3b-extract-compile-birdtree.sh` - bash script to iterate through all pieces of data and run `3b-extract-birdtree.sbatch`
    * `3b-extract-birdtree.sbatch` - bash script to load one piece of data and run `3b-extract-avg-within-range-birdtree.R` on HPCC
    * `3c-get-master-file-birdtree.R` - generate master file with rows for species and columns for environmental and life history data
    * `3c-get-master-file-birdtree.sbatch` - bash script to run `3c-get-master-file-birdtree.R`
    * `3d-env-var-range.R` - explore spatial variation in environmental variation across species ranges
  * `4-phylo-cor` - get phylogenetic correlation matrix
    * `4a-bird-get-con-tree-pieces.R` - calculate consensus trees for 10 chunks of 100 phylogenies for birds
    * `4b-bird-get-cor-matrix.R` - calculate final consensus tree and phylo cor matrix for birds 
  * `5-model/`
    * `5-gl-phylo-vint-berk-oe.R` - fit gl ~ env var
  * `6-rasterize` - rasterize ranges
    * `6a-create-rasters.R` - create .tif files birds
    * `6a-raster-comile.sh` - bash script to iterate through all pieces of data and run `6a-raster.sbatch`
    * `6a-raster.sbatch` - bash script to load one piece of data and run `6a-create-rasters.R` on HPCC
    * `6b-stack-rasters.R` - script to stack rasters and produce main .tif
  * `7-figures` - figures
  * `Model_files/` - Stan model files
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
    * `delta.csv` - archival $\delta_T$ and $\delta_P$ grid data 

## Contact Information

Casey Youngflesh - cyoungf@clemson.edu


