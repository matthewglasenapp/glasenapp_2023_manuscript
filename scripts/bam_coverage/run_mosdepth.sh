#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=mosdepth_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --time=168:00:00

module load python-3.6.5

mosdepth --by regions.bed --no-per-base -t 4 --fast-mode --thresholds 10,20,30,50,100 droe /hb/groups/pogson_group/dissertation/data/dedup_mapped_bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam

python3 plot-dist.py droe.mosdepth.global.dist.txt