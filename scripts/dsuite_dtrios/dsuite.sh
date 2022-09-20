#!/bin/bash
#SBATCH --job-name=d_suite
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=d_suite_%j.out
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=125GB
#SBATCH --time=168:00:00

# --JKnum sets the number of JackKnife blocks. A block size of 1Mb would require 922 blocks. 

Dsuite="/hb/groups/pogson_group/dissertation/software/Dsuite/Build/Dsuite"
Dsuite_path="/hb/groups/pogson_group/dissertation/software/"
filtered_vcf="/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls.g.vcf.gz"
SETS="/hb/home/mglasena/dissertation/data/dsuite/franciscanus/SETS.txt"

#$Dsuite Dtrios $filtered_vcf SETS.txt -t tree_lividus.nwk -o out --no-f4-ratio --JKnum 922

# Parallel option
DtriosParallel="/hb/groups/pogson_group/dissertation/software/Dsuite/utils/DtriosParallel"
cores=24

$DtriosParallel --dsuite-path $Dsuite_path --cores $cores -k 922 -t tree.nwk $SETS $filtered_vcf