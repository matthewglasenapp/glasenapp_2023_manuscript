#!/bin/bash
#SBATCH --job-name=single_locus_trees
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=single_locus_trees_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=24:00:00

fasta_alignment_dir=""

iqtree -S $fasta_alignment_dir --prefix loci