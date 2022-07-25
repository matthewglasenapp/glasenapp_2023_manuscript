#!/bin/bash
#SBATCH --job-name=iq_single_locus_parallel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=iq_single_locus_parallel_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=24:00:00

module load python-3.6.5
python3 iq_single_locus_parallel.py

