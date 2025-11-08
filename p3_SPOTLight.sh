#!/usr/bin/sh
#SBATCH --partition={partition_name}
#SBATCH -A {account_name}
#SBATCH --job-name=p3_SPOTLight
#SBATCH --time=120:00:00
#SBATCH -o Logs/p3_SPOTLight_%j.log
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1


echo Running on `hostname`
# source the required packages: R
export R_X=4.2.3

echo Started: `date`

Rscript p3.1_SPOTlight_refannolv3.r

echo Finished: `date`
