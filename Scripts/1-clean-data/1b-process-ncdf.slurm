#!/bin/bash
#SBATCH --job-name=ncdf
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --time=3:59:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kapsarke@msu.edu
#SBATCH --mem=400G #memory requested
#SBATCH --constraint=amd22
#SBATCH --output=%x-%j.SLURMout

#echos name of node
echo `hostname`

#load modules
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.2.0

#manually specify R lib
export R_LIBS=/mnt/ufs18/home-109/kapsarke/R/x86_64-pc-linux-gnu-library/4.2  

Rscript ./Scripts/1-clean-data/1b-process-ncdf.R /mnt/research/ibeem/variability_testing/data/L0/climate/era5/ /mnt/research/ibeem/variability_testing/data/L1/climate/era5/ 1,2,3,4,5,6,7,8,9,10,11,12
