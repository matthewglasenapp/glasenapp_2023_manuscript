#!/bin/bash

# Convert gff3 file to bed file using the convert2bed program of BEDOPS v2.4.39. 
# Instructions for download BEDOPS can be found here: https://bedops.readthedocs.io/en/latest/content/installation.html#installation

input_file=protein_coding_genes.gff
output_file=protein_coding_genes.bed

convert2bed --input=gff --output=bed < $input_file > $output_file
