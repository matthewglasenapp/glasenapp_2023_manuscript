#!/bin/bash
#SBATCH --job-name=phylonet
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=phylonet_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=125GB
#SBATCH --time=168:00:00

module load java/java-8
java -Xmx125g -jar PhyloNet_3.8.0.jar file.nex