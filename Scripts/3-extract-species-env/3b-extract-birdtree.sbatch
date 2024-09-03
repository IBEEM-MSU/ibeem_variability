#!/bin/bash
# Number of nodes needed:
#SBATCH --nodes=1
#
# Tasks per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=40G
#
# Wall time (e.g. "minutes", "hours:minutes:seconds", "days-hours", "days-hours:minutes"):
#SBATCH --time=3:59:00
#
# Mail type:
# #SBATCH --mail-type=FAIL,END
#
# Mail user:
# #SBATCH --mail-user=kapsarke@msu.edu
# #SBATCH --mail-user=ccy@msu.edu
# #SBATCH --account plzbuyin #Space lab node

# module load GCC/8.3.0
# module load OpenMPI/3.1.4
module load R-bundle-CRAN/2023.12-foss-2023a

#manually specify R lib
# export R_LIBS=/mnt/home/kapsarke/R/x86_64-pc-linux-gnu-library/4.1
# export R_LIBS=/mnt/home/ccy/R/x86_64-pc-linux-gnu-library/4.3

# Rscript ./Scripts/3-extract-species-env/3b-extract-avg-within-range.R $fileName 
Rscript ./Scripts/3-extract-species-env/3b-extract-avg-within-range-birdtree.R $fileName 

