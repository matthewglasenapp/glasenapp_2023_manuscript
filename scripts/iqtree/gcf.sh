#!/bin/bash
#SBATCH --job-name=gcf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=gcf_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=10GB
#SBATCH --time=24:00:00

species_tree_file=""
single_locus_trees_file=""
fasta_alignment_dir=""


iqtree -t $species_tree_file --gcf $single_locus_trees_file --cf-verbose --df-tree -p $fasta_alignment_dir --scf 100 --prefix concord -T 10