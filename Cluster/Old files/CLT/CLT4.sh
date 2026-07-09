#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=80:00:00
#SBATCH --job-name="CLT fits"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
R --slave CMD BATCH CLT4.R