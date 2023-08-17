#!/bin/bash

for i in {1..21}
do
  sbatch chunks/1c-chunk-$i.slurm
done
