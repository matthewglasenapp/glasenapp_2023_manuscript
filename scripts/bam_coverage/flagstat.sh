#!/bin/bash
#SBATCH --job-name=flagstat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=flagstat_%j.out
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

module load samtools/samtools-1.14
module load python-3.6.5

python3 -u flagstat.py