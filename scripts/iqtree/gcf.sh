#!/bin/bash
#SBATCH --job-name=gcf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=gcf_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0   
#SBATCH --time=7-0

species_tree_file="/hb/scratch/mglasena/strongylocentrotidae/species_tree/rev_dna.treefile"
single_locus_trees_file="/hb/scratch/mglasena/strongylocentrotidae/loci.treefile"
fasta_alignment_dir="/hb/scratch/mglasena/strongylocentrotidae/vcf2fasta_CDS/"

iqtree -t $species_tree_file --gcf $single_locus_trees_file --cf-verbose --df-tree -p $fasta_alignment_dir --scf 300 --prefix concord -T 8