#!/bin/sh
#SBATCH --account=stats
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=20
#SBATCH --time=400:00:00
#SBATCH --job-name="Cluster config 1"
#SBATCH --mail-user=Greg.Distiller@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
 
hostname
R --slave CMD BATCH Clusters1b.R