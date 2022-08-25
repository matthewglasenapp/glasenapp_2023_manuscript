#!/bin/bash
#SBATCH --job-name=d_suite
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=d_suite_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=24:00:00

# --JKnum sets the number of JackKnife blocks. A block size of 1Mb would require 922 blocks. 

Dtrios="/hb/groups/pogson_group/dissertation/software/Dsuite/BuildDsuite"
filtered_vcf="/hb/scratch/mglasena/data/combined_vcf/filtered_genotype_calls.g.vcf.gz"

$Dtrios $filtered_vcf SETS.txt -t tree.nwk -o out --no-f4-ratio --JKnum 922

# Parallel option
#DtriosParallel="/hb/groups/pogson_group/dissertation/software/Dsuite/utils/DtriosParallel"
#filtered_vcf="/hb/scratch/mglasena/data/combined_vcf/filtered_genotype_calls.g.vcf.gz"
#cores=24

#DtriosParallel --cores $cores -k 922 -t tree.nwk $filtered_vcf SETS.txt

