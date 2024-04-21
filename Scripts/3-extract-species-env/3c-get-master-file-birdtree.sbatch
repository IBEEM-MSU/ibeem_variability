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
#SBATCH --time=20:00:00
#
# Mail type:
#SBATCH --mail-type=FAIL,END
#
# Mail user:
# #SBATCH --mail-user=doserjef@msu.edu
#SBATCH --mail-user=kapsarke@msu.edu
# #SBATCH --account plzbuyin #Space lab node

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

#manually specify R lib
export R_LIBS=/mnt/home/kapsarke/R/x86_64-pc-linux-gnu-library/4.2

# Rscript /mnt/home/doserjef/ibeem_variability/Scripts/3-extract-species-env/3c-get-master-file-birdtree.R
Rscript ./Scripts/3-extract-species-env/3c-get-master-file-birdtree.R
