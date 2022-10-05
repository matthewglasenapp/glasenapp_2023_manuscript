#!/bin/bash
#SBATCH --job-name=phylonet
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=phylonet_%J.out
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0   
#SBATCH --time=20-0

module load java/java-8

phylonet="/hb/groups/pogson_group/dissertation/software/phylonet/PhyloNet_3.8.0.jar"

java -Xmx125g -jar $phylonet retic_0