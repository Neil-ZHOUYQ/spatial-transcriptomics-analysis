#!/usr/bin/sh
#SBATCH --partition={partition_name}
#SBATCH -A {account_name}
#SBATCH --job-name=p2_spatialtrans_integration
#SBATCH --time=120:00:00
#SBATCH -o Logs/p2_spatialtrans_integration_%j.log
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1

echo Running on `hostname`
# source the required packages: R
export R_X=4.2.3

echo Started: `date`

Rscript p2.1_spatialtrans_preprocessing.r
Rscript p2.2_spatialtrans_integration.r

echo Finished: `date`
