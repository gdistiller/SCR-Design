#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=25
#SBATCH --time=60:00:00
#SBATCH --job-name="GA4 Traps Both 40 (2nd set)"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.3.3-usr
srun R --slave CMD BATCH GABothStrataPars40_3c.R