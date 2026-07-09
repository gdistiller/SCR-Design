#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=5
#SBATCH --time=60:00:00
#SBATCH --job-name="Test Parallel"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
R --slave CMD BATCH TestParallel.R