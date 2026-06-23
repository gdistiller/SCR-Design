#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=25
#SBATCH --time=60:00:00
#SBATCH --job-name="GA Traps Both Max 40"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.3.3-usr
srun R --slave CMD BATCH GABothStrataParsMax40_3.R