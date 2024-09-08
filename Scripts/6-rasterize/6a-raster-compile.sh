#!/bin/bash
n="ras" # job name

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
#SBATCH --mem=5G
#
# Wall time (e.g. "minutes", "hours:minutes:seconds", "days-hours", "days-hours:minutes"):
#SBATCH --time=00:05:00

# Note that 98 is the number of pieces the BT IDs were split up into. 
# Should probably not hardcode this...
for i in {1..98}
do
  sbatch --job-name=$n.$i --output=$n.$i.SLURMout --export=fileName=BTIDsPiece-$i.rda ./Scripts/6-rasterize/6a-raster.sbatch
done
