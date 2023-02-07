# Get all coding exons from gff file and convert to bed file. Sort by exon name to keep hierarchical structure of gff3

wget -O - -o /dev/null https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_genomic.gff.gz | \
gunzip --stdout | \
awk '$3 == "exon"' | \
grep "gbkey=mRNA" | \
grep "ID=exon-NM\|ID=exon-XM" \
> exons.gff

convert2bed --input=gff --output=bed < exons.gff > exons.bed

cat exons.bed | sort -k4,4  > sorted_exons.bed