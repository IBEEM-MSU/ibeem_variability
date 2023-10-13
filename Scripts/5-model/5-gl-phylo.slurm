#!/bin/bash
#SBATCH --job-name=gl-phylo
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 6 #cpus
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ccy@msu.edu
#SBATCH --mem=80G #memory requested
#SBATCH --constraint=amd22
#SBATCH --output=%x-%j.SLURMout

#echos name of node
echo `hostname`

#load modules
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load R/4.1.0

#manually specify R lib
export R_LIBS=/mnt/home/ccy/R/x86_64-pc-linux-gnu-library/4.1

Rscript /mnt/home/ccy/variability/Scripts/5-model/5-gl-phylo.R
