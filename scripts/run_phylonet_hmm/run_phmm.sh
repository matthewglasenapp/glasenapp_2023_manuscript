#!/bin/bash
#SBATCH --job-name=phylonet_hmm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=phmm_%J.out
#SBATCH --error=phmm_%J.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=15GB
#SBATCH --array=[0-20]
#SBATCH --time=7-0

module load python-3.6.5

array_id=$SLURM_ARRAY_TASK_ID
export array_id

python3 -u /hb/home/mglasena/dissertation/scripts/phylonet_hmm/run_phylonet_hmm.py
