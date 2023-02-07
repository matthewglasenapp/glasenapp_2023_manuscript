The scrpt preprocess.py pre-processes a list of raw fastq file pairs. Adapter sequences are marked, reads are aligned to the S. purpuratus reference genome, duplicate sequences are marked, and variants are called. The script preprocess.sh runs preprocess.py as a batch array job on UCSC Hummingbird. 

The script genotype_vcfs.py merges single vcf files to create a multisample vcf file and performs joint genotyping. 

The script filter_genotypes.py filters the raw genotype calls. 