# IBEEM_variability

## Project Links:

[Project Google Drive](https://drive.google.com/drive/u/1/folders/1hnJP1CRZQSph2M2cIhCEwvujxFOhfM0k)

[Meeting Notes](https://docs.google.com/document/d/1ix5mSCbO7rRCj2juQN04UeabejayMgcYT0-yo08uSZ0/edit)

[Data Documentation](https://docs.google.com/document/d/13Hn0klwabOQzdfCg1W-FIF00j1t_YgT40v-99Fhkt3Y/edit)

[Proposal](https://drive.google.com/file/d/1K0jelDSM3ZlyHDlNI3SRYAbzlrS0hpiU/view?usp=share_link)


## Repo structure:

* `Scripts/`
  * `0-clean-Bird-et-al.R` - clean and combine gen time data from Bird et al. 2020
  * `1-process-ncdf.R` - convert ERA netcdf temp and precip to csvs
  * `2-env-metrics.R` - read in env data, calculate env metrics, and plot
  * `example_netcdf_processing.R` - netcdf processing examples
* `Sample_output/` - small data objects from env variable processing (to work with data and avoid reprocessing)
* `Data/` (ignored)

## Calc Env var:

* Modules to be loaded for interactive session (coudl batch as well):
  * `module load GCC/11.2.0`
  * `module load OpenMPI/4.1.1`
  * `module load R/4.1.2`
* Request high mem interactive session
  * `salloc -N 1 -c 4 --time=3:59:00 --constraint=amd22 --mem=400gb`
* Average ER5 data over specified months (year, and JJA) to produce L1 data (3 args: in dir, out dir, months)
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 1,2,3,4,5,6,7,8,9,10,11,12`
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/1-process-ncdf.R /mnt/research/ibeem/L0/climate/era5/ /mnt/research/ibeem/L1/climate/era5/ 6,7,8`
* Calc env variability metrics to produce L2 data (2 args: in dir, out dir)
  * `Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics.R /mnt/research/ibeem/L1/climate/era5/ /mnt/research/ibeem/L2/climate/era5/`

## Data notes:
* Env L0 - raw ERA data (`/mnt/research/ibeem/L0/climate/era5`)
* Env L1 - ERA data averaged to specified months (`/mnt/research/ibeem/L1/climate/era5`)
* Env L2 - env variability metrics (`/mnt/research/ibeem/L2/climate/era5`)
