#!/bin/bash
#SBATCH --job-name=create_phmm_input
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=create_phmm_input_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=2-0

module load python-3.6.5

python3 -u /hb/home/mglasena/dissertation/scripts/phylonet_hmm/scaffolds.py
python3 -u /hb/home/mglasena/dissertation/scripts/phylonet_hmm/create_matrices.py