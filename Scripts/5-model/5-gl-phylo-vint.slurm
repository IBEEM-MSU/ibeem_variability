#!/bin/bash
#SBATCH --job-name=gl-pv-bird
#SBATCH -N 1 #number of tasks
#SBATCH -n 1 #number of nodes
#SBATCH -c 6 #cpus
#SBATCH --time=168:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kapsarke@msu.edu
#SBATCH --mem=80G #memory requested
#SBATCH --constraint=amd22
#SBATCH --output=%x-%j.SLURMout

#echos name of node
echo `hostname`

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load R/4.1.2

#manually specify R lib
export R_LIBS=/mnt/home/kapsarke/R/x86_64-pc-linux-gnu-library/4.1

Rscript ./Scripts/5-model/5-gl-phylo-vint.R
