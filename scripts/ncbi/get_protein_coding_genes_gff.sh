# Get all protein coding genes from gff file in gff3 format. 
wget -O - -o /dev/null https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_genomic.gff.gz | \
gunzip --stdout | \
awk '$3 == "gene"' | \
grep -i "gene_biotype=protein_coding" \
> protein_coding_genes.gff