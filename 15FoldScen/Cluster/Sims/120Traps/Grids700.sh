#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=40
#SBATCH --time=400:00:00
#SBATCH --job-name="Grids 700 Sims 120"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.5.2
R --slave CMD BATCH Grids700.R