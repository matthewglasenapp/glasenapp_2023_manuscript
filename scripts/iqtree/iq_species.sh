#!/bin/bash
#SBATCH --job-name=iq_species
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=iq_species.out
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0   
#SBATCH --time=7-0

fasta_alignment_dir="/hb/scratch/mglasena/strongylocentrotidae/vcf2fasta_CDS/"
iqtree="/hb/groups/pogson_group/dissertation/software/iqtree-2.2.2.2/bin/iqtree2"

iqtree -p $fasta_alignment_dir -m MFP -B 1000 --sampling GENESITE -T 8 --prefix rev_dna