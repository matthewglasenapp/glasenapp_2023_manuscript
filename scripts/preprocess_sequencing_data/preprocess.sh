#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=preprocess_%J.out
#SBATCH --error=preprocess_%J.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60GB
#SBATCH --array=[0-9]%6

module load python-3.6.5
module load samtools/samtools-1.14
module load bcftools/bcftools-1.14

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u preprocess.py