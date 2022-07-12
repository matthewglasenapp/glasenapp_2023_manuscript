#!/bin/bash
#SBATCH --job-name=fasterq_dump
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=fasterq_dump_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=168:00:00

module load python-3.6.5
python3 -u /hb/groups/pogson_group/dissertation/scripts/download_sra_data/fasterq_dump.py