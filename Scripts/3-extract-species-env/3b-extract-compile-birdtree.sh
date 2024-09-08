#!/bin/bash
n="es" # job name

# Note that 98 is the number of pieces the BT IDs were split up into. 
# Should probably not hardcode this...
for i in {1..98}
do
  sbatch --job-name=$n.$i --output=$n.$i.SLURMout --export=fileName=BTIDsPiece-$i.rda ./Scripts/3-extract-species-env/3b-extract-birdtree.sbatch
done
