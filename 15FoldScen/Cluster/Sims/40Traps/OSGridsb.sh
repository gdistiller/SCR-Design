#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=20
#SBATCH --time=400:00:00
#SBATCH --job-name="Optimal Grid"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.3.3-usr
R --slave CMD BATCH OSGridsb.R