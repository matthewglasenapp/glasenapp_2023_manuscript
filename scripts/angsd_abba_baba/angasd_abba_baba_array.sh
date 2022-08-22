#!/bin/bash
#SBATCH --job-name=angsd_abba_baba
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=angsd_%J.out
#SBATCH --error=angsd_%J.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=[0-5]

module load python-3.6.5

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u angsd_abba_baba.py