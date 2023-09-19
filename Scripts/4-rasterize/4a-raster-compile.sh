#!/bin/bash
n="ras" # job name

# Note that 110 is the number of pieces the BL IDs were split up into. 
# Should probably not hardcode this...
for i in {1..110}
do
  sbatch --job-name=$n.$i --output=$n.$i.SLURMout --export=fileName=BLIDsPiece-$i.rda 4a-raster.sbatch
done
