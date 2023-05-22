#!/bin/bash
#SBATCH --job-name=env-GAM-6_8
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --time=60:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ccy@msu.edu
#SBATCH --mem=50G #memory requested
#SBATCH --output=%x-%j.SLURMout
#SBATCH --account plzbuyin #Space lab node

#echos name of node
echo `hostname`

#load modules
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.1.0

#manually specify R lib
export R_LIBS=/mnt/home/ccy/R/x86_64-pc-linux-gnu-library/4.1

#months 6-8 GAM
Rscript /mnt/research/ibeem/ibeem_variability/Scripts/2-env-metrics/2-env-metrics-GAM.R /mnt/research/ibeem/L1/climate/era5/ERA5-6_7_8.csv /mnt/research/ibeem/L2/climate/era5/Env-var-GAM-6_7_8.csv