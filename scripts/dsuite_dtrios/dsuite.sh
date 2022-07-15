#!/bin/bash
#SBATCH --job-name=d_suite
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=d_suite_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=24:00:00

filtered_vcf="/hb/groups/pogson_group/dissertation/data/combined_vcf/genotype_calls.g.vcf.gz"
/hb/home/mglasena/software/Dsuite/Build/Dsuite Dtrios $filtered_vcf SETS.txt -t tree.txt -o out --no-f4-ratio