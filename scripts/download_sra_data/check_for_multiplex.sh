#!/bin/bash
#SBATCH --job-name=check_for_multiplex.py
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=check_for_multiplex_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=10GB
#SBATCH --time=24:00:00

module load python-3.6.5
python3 -u check_for_multiplex.py