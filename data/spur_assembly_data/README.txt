This directory contains some manipulated files from the S. purpuratus reference assembly. 

The file unique_exons.bed is a bed file of all unique, protein-coding exons in the S. purpuratus assembly. Duplicate exons (exons with different names but identical start and stop coordinates) were removed.

The file sorted_exons.bed is a bed file of all protein-coding exons in the S. purpuratus assembly. The bed file is sorted by exon name to keep the hierarchical structure of gff3. There are duplicate exons in cases where an exon is part of more than one transcript. Duplicate exons were retained because the coverage metrics for these exons were aggregated to calculate coverage metrics for each mRNA record.

The file protein_coding_genes.bed is a bed file of all protein-coding gene records in the S. purpuratus assembly. 

The file protein_coding_genes.gff is a ggf3 file of all protein-coding gene records in the S. purpuratus assembly. 