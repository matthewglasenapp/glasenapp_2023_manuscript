#!/bin/bash
#SBATCH --job-name=iq_species
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=iq_species_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=24:00:00

fasta_alignment_dir=""

iqtree -p $fasta_alignment_dir -B 1000 -T 10 --prefix rev_dna