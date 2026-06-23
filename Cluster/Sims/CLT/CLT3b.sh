#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=2 --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --job-name="CLT 1 Occ (comp mask)"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
R --slave CMD BATCH CLT3b.R