#!/bin/bash
#SBATCH --job-name=angsd_abbababa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=andsd_abbababa_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --mem=60GB
#SBATCH --time=72:00:00

module load R/R-3.6.1

reference="/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

angsd -doAbbababa 1 \
-doCounts 1 \
-baq 1 \
-ref $reference \
-useLast 1 \
-bam bam.filelist \
-out out \
-blockSize 1000000 \
-minMapQ 20 \
-minQ 20 \
-only_proper_pairs 1 \
-remove_bads 1 \
-uniqueOnly 1 \


Rscript jackKnife.R file=out.abbababa indNames=bam.filelistnames outfile=results

######
#Command line options copy/pasted from ANGSD documentation

#-doAbbababa	Perform an ABBA-BABA test
#-doCounts	Calculate various counts statistics

#baq adjust qscores around indels

#-useLast [int]
#1: use the last individual in the bam file list as outgroup instead of a fasta file (-anc)

#-enhance [int]
#1: use only sites where the reads for the outgroup has the same base for all reads. Only works with -useLast 1

#-minMapQ [int]=0
#Minimum mapQ quality.

#-minQ [int]=0
#Minimum base quality score.

#-only_proper_pairs [int]=1
#Include only proper pairs (pairs of read with both mates mapped correctly). 1: include only proper (default), 0: use all reads. Only relevant for paired end data.

#-remove_bads [int]=1
#Same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads). 0 no , 1 remove (default).

#-uniqueOnly [int]=0
#Remove reads that have multiple best hits. 0 no (default), 1 remove.

#-setMaxDepthInd	-1	(If depth persample is larger then individual is removed from analysis (from site).
#-setMinDepthInd	-1	(If depth persample is smaller then individual is removed from analysis (from site).

#Use -nThreads or -P for number of threads allocated to the program