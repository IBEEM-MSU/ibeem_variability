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
#SBATCH --mem=25G
#
# Wall time (e.g. "minutes", "hours:minutes:seconds", "days-hours", "days-hours:minutes"):
#SBATCH --time=10:00:00
#
# Mail type:
#SBATCH --mail-type=FAIL,END
#
# Mail user:
# #SBATCH --mail-user=kapsarke@msu.edu
#SBATCH --mail-user=ccy@msu.edu
#SBATCH --account plzbuyin #Space lab node

module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.0.2

Rscript /mnt/home/ccy/variability/Scripts/3-extract-species-env/3b-extract-avg-within-range-terrestrial-mammal.R $fileName 

