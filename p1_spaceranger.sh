#!/usr/bin/sh
#SBATCH --partition={partition_name}
#SBATCH -A {account_name}
#SBATCH --job-name=p1_spanceranger
#SBATCH --time=120:00:00
#SBATCH -o Logs/p1_spanceranger_%j.log
#SBATCH --mem=150G
#SBATCH --cpus-per-task=16
#SBATCH --time=120:00:00

echo Running on `hostname`
# source the required packages: spaceranger
echo Started: `date`

reference=/path_to_reference/
sample_folder=/path_to_ST_sample/
image_file=/path_to_ST_sample_image/
probeset_file=/path_to_probeset/

# Run spaceranger count
spaceranger count --id="id_name" \
                  --description="A short description" \
                  --transcriptome=$reference \
                  --fastqs=$sample_folder \
                  --cytaimage=$image_file \
                  --probe-set=$probeset_file \
                  --slide=V43L10-047 \
                  --area=D1 \
                  --localcores=16 \
                  --localmem=128

echo "Finished"
