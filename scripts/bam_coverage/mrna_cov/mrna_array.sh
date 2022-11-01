#!/bin/bash
#SBATCH --job-name=mrna_cov
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=mrna_cov_%J.out
#SBATCH --error=mrna_cov_%J.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=[0-12]
#SBATCH --time=20-0 

module load python-3.6.5
module load R

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u mrna.py