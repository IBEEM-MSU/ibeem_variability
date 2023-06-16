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
    * `1b-process-ncdf.R` - create time series of env data averaged over specified months
    * `1c-taxonomic-harmonization.R` - sort out naming differences
  * `2-env-metrics/`
    * `2-env-metrics.R` - calculate env variability metrics
    * `2-env-metrics-GAM.R` - calculate env variability metrics using GAM detrend
  * `3-explore-env-var.R` - explore environmental variability metrics
  * `example_netcdf_processing.R` - netcdf processing examples
* `Sample_output/` - small data objects from env variable processing (to work with data and avoid reprocessing)
* `Data/` (ignored)

## Data sources:
* [ERA5 reanalysis data (temp and precipitation)](https://rda.ucar.edu/datasets/ds633.1/)
* [Bird Life International range maps](http://datazone.birdlife.org/species/requestdis)
* [IUCN mammal range amps](https://www.iucnredlist.org/resources/spatial-data-download)
* [Bird et al. 2020 *Conservation Biology* bird generation time](https://conbio.onlinelibrary.wiley.com/doi/10.1111/cobi.13486)
* [Pacifici et al. 2013 *Nature Conservation* mammal generation time](https://natureconservation.pensoft.net/article/1343/download/pdf/)
* [Elton Traits (Wilman et al. 2014) for birds and mammals](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/13-1917.1)
* [PanTHERIA (Jones et al. 2009) for mammals](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/08-1494.1)
* [AVONET (Tobias et al. 2022) for birds](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898)
* [Coonety et al. 2020 for birds](https://www.nature.com/articles/s41467-020-16257-x)

## Data processing levels:
* Env L0 - raw ERA data (`/mnt/research/ibeem/L0/climate/era5/`)
* Env L1 - ERA data averaged to specified months (`/mnt/research/ibeem/L1/climate/era5/`)
* Env L2 - env variability metrics (`/mnt/research/ibeem/L2/climate/era5/`)

## Env time series (avg over specified months):
* Request high mem interactive session (could batch as well) - may not need nearly this much memory:
  * `salloc -N 1 -c 4 --time=3:59:00 --constraint=amd22 --mem=400gb`
* Modules to be loaded:
  * `module load GCC/8.3.0`
  * `module load OpenMPI/3.1.4`
  * `module load R/4.1.0`
* Average ER5 data over specified months (year, and JJA) to produce L1 data (3 args: in dir, out dir, months):
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-clean-data/1b-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 1,2,3,4,5,6,7,8,9,10,11,12`
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-clean-data/1b-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 6,7,8`

## Env variability metrics:
* Request interactive session (could/should batch - see `Scripts/2-env-metrics/env-1_12.slurm` and `Scripts/2-env-metrics/env-6_8.slurm`):
  * `salloc -N 1 -c 4 --time=40:00:00 --mem=50gb`
* Modules to be loaded:
  * `module load GCC/8.3.0`
  * `module load OpenMPI/3.1.4`
  * `module load R/4.1.0`
* Calc env variability metrics to produce L2 data (2 args: in file, out file):
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics/2-env-metrics.R /mnt/research/ibeem/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv /mnt/research/ibeem/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv`
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics/2-env-metrics.R /mnt/research/ibeem/L1/climate/era5/ERA5-6_7_8.csv /mnt/research/ibeem/L2/climate/era5/Env-var-6_7_8.csv`

## Extracting env and trait data for each species
+ Code for initial processing for birds (breeding season) is in `Scripts/4-extract-species-env`
+ Currently using species with ITIS ids, did not resolve NAs. 
+ Generates a main file with each row a species (with associated ITIS ID and accepted name). Each column contains either the average of one of the environmental variables across the range of the species, or some trait value. 
+ Main file is stored in `/mnt/research/ibeem/data/L2/main-bird-data.csv` and is also in the [`L2` folder on Google drive](https://drive.google.com/drive/folders/1c4dF8AEgOf7zvVUd5GTUT7vzj5hoLuGj). 
+ Trait data in the master file currently include values from Bird et al and AVONET. Can add in ELTON later on if there is something specific we want from there. The column name has the data source as a prefix. See metadata for the data sources for specific info.

