#!/bin/bash
#SBATCH --job-name=ncdf-monthly-4
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 2 #cpus
#SBATCH --time=3:59:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ccy@msu.edu
#SBATCH --mem=50G #memory requested
#SBATCH --output=%x-%j.SLURMout
#SBATCH --constraint=amd22

#echos name of node
echo `hostname`

#load modules
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.1.0

#manually specify R lib
export R_LIBS=/mnt/home/ccy/R/x86_64-pc-linux-gnu-library/4.1

Rscript /mnt/research/ibeem/variability/Scripts/1-clean-data/1c-process-ncdf-monthly.R /mnt/research/ibeem/variability/data/L0/climate/era5/ /mnt/research/ibeem/variability/data/L1/climate/era5/ 4

