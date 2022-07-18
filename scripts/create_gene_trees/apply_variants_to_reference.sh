#!/bin/bash
#SBATCH --job-name=array_consensus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=consensus_%J.out
#SBATCH --error=consensus_%J.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[0-8]

module load python-3.6.5
module load bcftools/bcftools-1.14

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u apply_variants_to_reference.py