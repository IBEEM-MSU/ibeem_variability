# Variability_life_history

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species. 

## Table of Contents
- [Overview](#Overview)
- [Funding](#Funding)
- [Associated publications](#Associated-publications)
- [Data](#Data)
- [Spatiotemporal extent and resolution](#Spatiotemporal-extent-and-resolution)
- [Workflow](#Workflow)
- [Contact Information](#Contact-information)

## Overview

This repository contains code to assess the relationship between environmental variability and pace of life across the world's resident bird species. The project is part of the [Institute for Biodiversity, Ecology, Evolution, and Macrosystems (IBEEM)](https://ibeem.msu.edu/) at Michigan State University, and led by IBEEM Postdoctoral Fellows. The research uses data from open access biodiversity repositories to evaluate how patterns of environmental variability within a year (i.e., seasons) and between years affect species' life histories. The analysis investigates environmental variability within and between years across the ranges of 7,477 non-migratory and non-marine avian species in order to evaluate the impacts of temporal environmental variability on species' pace of life. Environmental data are from the ERA5 climate reanalysis, bird trait and distribution data are from AVONET, BirdTree, and BirdLife databases as well as previously published research. Environmental variability values were extracted across individual species' ranges in order to quantify how environmental change has affected avian species' pace of life around the world. For more on the project see the [project research summary](https://ibeem.msu.edu/research--news.html) and the publications linked below. 

## Funding

Funding is provided by a Michigan State University Strategic Partnership Grant for the [Institute for Biodiversity, Ecology, Evolution, and Macrosystems (IBEEM)](https://ibeem.msu.edu/).

## Associated publications

Youngflesh, C, K Kapsar, A Uscanga, PJ Williams, JW Doser, L Kounta, PL Zarnetske. Environmental variation shapes life history of the world’s birds. *In Review*

## Data

* Youngflesh, C., K. Kapsar, A. Uscanga, P.J. Williams, J.W. Doser, L. Kounta, and P.L. Zarnetske. 2024. Inter- and intra-annual temperature and precipitation variability (1950-2022) across the ranges of non-migratory birds and their association with generation length ver 2. Environmental Data Initiative. [DOI: https://doi.org/10.6073/pasta/fc057e288924cb09fc3eb5092856d99c](https://doi.org/10.6073/pasta/fc057e288924cb09fc3eb5092856d99c) (Accessed 2024-09-11).
  * `env_var.csv` - Global grid of environmental variability
  * `delta.csv` - Global grid of environmental change scaled by variability and generation length 
  * `final-bird-data-for-archival.csv` - Model data: Species generation length as a function of environmental variability
  * `bird-consensus-tree.rda` - Consensus phylogenetic tree

* Bird, J. P., R. Martin, H. R. Akçakaya, J. Gilroy, I. J. Burfield, S. Garnett, A. Symes, J. Taylor, Ç. H. Şekercioğlu, and S. H. M. Butchart. 2020. Generation lengths of the world’s birds and their implications for extinction risk. Conservation Biology. 34(5), 1252-1261. [DOI: doi.org/10.1111/cobi.13486](doi.org/10.1111/cobi.13486)
  * `cobi13486-sup-0003-tables3.xlsx` (Supplemental table)
  * `cobi13486-sup-0004-tables4.xlsx` (Supplemental table)

* Tobias, J. A., C. Sheard, A. L. Pigot, A. J. M. Devenish, J. Yang, F. Sayol, M. H. C. Neate‐Clegg, N. Alioravainen, T. L. Weeks, R. A. Barber, P. A. Walkden, H. E. A. MacGregor, S. E. I. Jones, C. Vincent, A. G. Phillips, N. M. Marples, F. A. Montaño‐Centellas, V. Leandro‐Silva, S. Claramunt, B. Darski, B. G. Freeman, T. P. Bregman, C. R. Cooney, E. C. Hughes, E. J. R. Capp, Z. K. Varley, N. R. Friedman, H. Korntheuer, A. Corrales‐Vargas, C. H. Trisos, B. C. Weeks, D. M. Hanz, T. Töpfer, G. A. Bravo, V. Remeš, L. Nowak, L. S. Carneiro, A. J. Moncada R., B. Matysioková, D. T. Baldassarre, A. Martínez‐Salinas, J. D. Wolfe, P. M. Chapman, B. G. Daly, M. C. Sorensen, A. Neu, M. A. Ford, R. J. Mayhew, L. Fabio Silveira, D. J. Kelly, N. N. D. Annorbah, H. S. Pollock, A. M. Grabowska‐Zhang, J. P. McEntee, J. Carlos T. Gonzalez, C. G. Meneses, M. C. Muñoz, L. L. Powell, G. A. Jamie, T. J. Matthews, O. Johnson, G. R. R. Brito, K. Zyskowski, R. Crates, M. G. Harvey, M. Jurado Zevallos, P. A. Hosner, T. Bradfer‐Lawrence, J. M. Maley, F. G. Stiles, H. S. Lima, K. L. Provost, M. Chibesa, M. Mashao, J. T. Howard, E. Mlamba, M. A. H. Chua, B. Li, M. I. Gómez, N. C. García, M. Päckert, J. Fuchs, J. R. Ali, E. P. Derryberry, M. L. Carlson, R. C. Urriza, K. E. Brzeski, D. M. Prawiradilaga, M. J. Rayner, E. T. Miller, R. C. K. Bowie, R. Lafontaine, R. P. Scofield, Y. Lou, L. Somarathna, D. Lepage, M. Illif, E. L. Neuschulz, M. Templin, D. M. Dehling, J. C. Cooper, O. S. G. Pauwels, K. Analuddin, J. Fjeldså, N. Seddon, P. R. Sweet, F. A. J. DeClerck, L. N. Naka, J. D. Brawn, A. Aleixo, K. Böhning‐Gaese, C. Rahbek, S. A. Fritz, G. H. Thomas, and M. Schleuning. 2022. AVONET: morphological, ecological and geographical data for all birds. Ecology Letters 25:581–597. [DOI: https://doi.org/10.1111/ele.13898](https://doi.org/10.1111/ele.13898)
  * `AVONET1_Birdlife.csv` (Tab 2 of `Supplementary dataset 1.xlsx`)
  * `avonet-birdtree-crosswalk.csv` (Tab 11 of `Supplementary dataset 1.xlsx`)

* European Centre for Medium-Range Weather Forecasts. 2022. ERA5 monthly mean back extension 1950-1978 (Preliminary version). Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory. (DOI: [https://doi.org/10.5065/JAXB-X906]([https://doi.org/10.5065/JAXB-X906))
  * `.nc` files for monthly temperature and precipitation

* European Centre for Medium-Range Weather Forecasts. 2019, updated yearly. ERA5 Reanalysis (Monthly Mean 0.25 Degree Latitude-Longitude Grid). Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory. (DOI: [https://doi.org/10.5065/P8GT-0R61](https://doi.org/10.5065/P8GT-0R61))
  * `.nc` files for monthly temperature and precipitation

* BirdLife International and Handbook of the Birds of the World. (2022). Bird species distribution maps of the world. Version 2022.2. Available at [http://datazone.birdlife.org/species/requestdis](http://datazone.birdlife.org/species/requestdis).
  * `BOTW.gdb/a0000000a.gdbtable`
  
* Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. Nature 491:444–448.(DOI: [https://doi.org/10.1038/nature11631](https://doi.org/10.1038/nature11631))
  * `.tre` files from `HackettStage2_XXX.zip` downloads

## Spatiotemporal extent and resolution

- Spatial extent: global
- Spatial resolution: 0.25° (climate data)
- Temporal extent: 1950–2022 (climate data)

## Usage

Software used to analyze data includes: R version 4.3.2

## Workflow

The workflow for this repository follows the guidelines set out by the Environmental Data Initiative (EDI). Briefly, this involves aligning with FAIR data practices, and employing a workflow that uses different levels for harmonization and derived data products. Data are read in as raw data at Level 0 (L0). Level 1 (L1) data represent cleaned L0 data. Level 2 data are data merged or otherwise derived from two or more L1 data, etc. Raw data citations with associated DOIs are given above. Data must be organized according to the below Repository structure.

Scripts are designed to be run in sequence from `0-XXXX.R` -> `1-XXXX.R` -> ... The directory where the data are located (i.e., the parent directory of `data/`) as well as the run date of `5-gl-phylo-vint-berk-oe.R` (for scripts `6-XXXX` and beyond) need to be specified in `0-config.R` before running. In some cases, High Performance Computing (HPC) resources are needed to processes data due to substantial compute requirements. This is the case for any script that has a `.slurm` script associated with it (e.g., `1b-process-ncdf.R` and `1b-process-ncdf.slurm`). `.slurm` scripts should be changed according to the specifics of the user's HPC system and the location of the data.



**Repository structure:**

* `Scripts/`
  * `0-config.R` - configuration file specifying location of data directory and model run date. `data/` parent directory and run date for `5-gl-phylo-vint-berk-oe.R` must be specified in this file as it will be unique to each user.
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
* `data/` (ignored)
  * `L0/` - raw data
    * `Bird_et_al_2020/`
      * `cobi13486-sup-0003-tables3.xlsx` (source: Bird et al. 2020)
      * `cobi13486-sup-0004-tables4.xlsx` (source: Bird et al. 2020)
    * `climate/era5/` - raw ERA5 reanalysis data
      * `T2m/` - temperature (.nc files; source for 1950-1978: https://doi.org/10.5065/JAXB-X906, 1979-2022: https://doi.org/10.5065/P8GT-0R61)
      * `PREC/` - precipitation (.nc files; source for 1950-1978: https://doi.org/10.5065/JAXB-X906, 1979-2022: https://doi.org/10.5065/P8GT-0R61)
    * `ranges/BOTW.gdb` - BirdLife range maps (source: BirdLife International and Handbook of the Birds of the World. 2022)
    * `trait/` - raw trait data
      * `Bird_et_al_gen_length_birds.csv` - processed from `1a-clean-Bird-et-al.R`
    * `phylogeny/BirdzillaHackett6.tre` - bird phylogeny (source: Jetz et al. 2012)
  * `L1/`
    * `climate/era5/` - ERA data averaged over specified months (one value per cell/year)
      * `ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv` - yearly average over all months and sd across months (seasonality) for temp and (sqrt root transform of) precip
    * `range/bird-breeding/` - bird ranges
      * `.shp` (and associated) files with breeding range for each species
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


