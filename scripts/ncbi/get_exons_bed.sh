# Get unique exons from gff file and convert to bed file

wget -O - -o /dev/null https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_genomic.gff.gz | \
gunzip --stdout | \
awk '$3 == "exon"' | \
grep "gbkey=mRNA" | \
grep "ID=exon-NM\|ID=exon-XM" | \
sort -u -k1,1 -k4,4n -k5,5n > unique_exons.gff

convert2bed --input=gff --output=bed < unique_exons.gff > unique_exons.bed