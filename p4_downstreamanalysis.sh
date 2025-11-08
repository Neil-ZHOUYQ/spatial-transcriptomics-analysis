#!/usr/bin/sh
#SBATCH --partition={partition_name}
#SBATCH -A {account_name}
#SBATCH --job-name=p3_downstreamanalysis
#SBATCH --time=120:00:00
#SBATCH -o Logs/p3_downstreamanalysis_%j.log
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1


echo Running on `hostname`
# source the required packages: R
export R_X=4.2.3

echo Started: `date`

Rscript p4.1_downstreamanalysis_degs.r

echo Finished: `date`
