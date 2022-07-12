#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=demultiplex_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=30GB
#SBATCH --time=168:00:00

module load python-3.6.5
python3 -u demultiplex.py