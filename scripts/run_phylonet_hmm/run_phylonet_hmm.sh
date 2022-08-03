#!/bin/bash
#SBATCH --job-name=phylonet_hmm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=phylonet_hmm_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=48:00:00

module load python-3.6.5

python3 /hb/home/mglasena/dissertation/scripts/phylonet/scaffolds.py
python3 /hb/home/mglasena/dissertation/scripts/phylonet/create_matrices.py
python3 /hb/home/mglasena/dissertation/scripts/phylonet/run_phylonet_hmm.py
python3 /hb/home/mglasena/dissertation/scripts/phylonet/process_hmm.py