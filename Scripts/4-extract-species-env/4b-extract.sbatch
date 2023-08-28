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
#SBATCH --mail-type=FAIL,END
#
# Mail user:
# #SBATCH --mail-user=doserjef@msu.edu
#SBATCH --mail-user=ccy@msu.edu
# #SBATCH --account plzbuyin #Space lab node

module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.1.0

Rscript /mnt/home/ccy/variability/Scripts/4-extract-species-env/4b-extract-avg-within-range.R $fileName 

