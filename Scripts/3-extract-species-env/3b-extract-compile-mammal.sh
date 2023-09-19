#!/bin/bash
n="es" # job name

# Note that 54 is the number of pieces the MAM IDs were split up into. 
# Should probably not hardcode this...
for i in {31..54}
do
  sbatch --job-name=$n.$i --output=$n.$i.SLURMout --export=fileName=MAMIDsPiece-$i.rda 3b-extract-mammal.sbatch
done
