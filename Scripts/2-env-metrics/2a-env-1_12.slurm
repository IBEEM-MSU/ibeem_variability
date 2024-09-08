#!/bin/bash
#SBATCH --job-name=env-1_12
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 4 #cpus
#SBATCH --time=48:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kapsarke@msu.edu
#SBATCH --mem=60G #memory requested
#SBATCH --output=%x-%j.SLURMout

#echos name of node
echo `hostname`

#load modules
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.2.0

#manually specify R lib
export R_LIBS=/mnt/ufs18/home-109/kapsarke/R/x86_64-pc-linux-gnu-library/4.2  

#months 1-12
Rscript ./Scripts/2-env-metrics/2a-env-metrics.R /mnt/research/ibeem/variability_testing/data/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv /mnt/research/ibeem/variability_testing/data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv
