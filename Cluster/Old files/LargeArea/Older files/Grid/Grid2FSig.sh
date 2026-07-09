#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=10
#SBATCH --time=60:00:00
#SBATCH --job-name="Sims Grid 2FSigma"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.3.3-usr
R --slave CMD BATCH Grid2FSig.R