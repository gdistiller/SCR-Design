#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=40
#SBATCH --time=400:00:00
#SBATCH --job-name="GA4 Traps Avg 120"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
module load R/R-4.3.3-usr
srun R --slave CMD BATCH GA4AvgStrataPars120_1b.R