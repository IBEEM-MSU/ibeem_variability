# IBEEM_variability

## Project Links:

[Project Google Drive](https://drive.google.com/drive/u/1/folders/1hnJP1CRZQSph2M2cIhCEwvujxFOhfM0k)

[Meeting Notes](https://docs.google.com/document/d/1ix5mSCbO7rRCj2juQN04UeabejayMgcYT0-yo08uSZ0/edit)

[Data Documentation](https://docs.google.com/document/d/13Hn0klwabOQzdfCg1W-FIF00j1t_YgT40v-99Fhkt3Y/edit)

[Proposal](https://drive.google.com/file/d/1K0jelDSM3ZlyHDlNI3SRYAbzlrS0hpiU/view?usp=share_link)


## Repo structure:

* `Scripts/`
  * `0-clean-Bird-et-al.R` - clean and combine gen time data from Bird et al. 2020
  * `1-process-ncdf.R` - create time series of env data averaged over specified months
  * `2-env-metrics.R` - calculate env variability metrics
  * `3-taxonomic-harmonization.R` - sort out naming differences
  * `X-explore-env-var.R` - explore environmental variabilty metrics
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

## Calculating Env var:
* Request high mem interactive session (could batch as well):
  * `salloc -N 1 -c 4 --time=3:59:00 --constraint=amd22 --mem=400gb`
* Modules to be loaded:
  * `module load GCC/11.2.0`
  * `module load OpenMPI/4.1.1`
  * `module load R/4.1.2`
* Average ER5 data over specified months (year, and JJA) to produce L1 data (3 args: in dir, out dir, months):
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 1,2,3,4,5,6,7,8,9,10,11,12`
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 6,7,8`
* Calc env variability metrics to produce L2 data (2 args: in file, out file):
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics.R /mnt/research/ibeem/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv /mnt/research/ibeem/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv`
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics.R /mnt/research/ibeem/L1/climate/era5/ERA5-6_7_8.csv /mnt/research/ibeem/L2/climate/era5/Env-var-6_7_8.csv`
