# Create conda environment on UCSC hummingbird for genomic analyses

module load miniconda3.9

conda create -y -n analysis

conda activate analysis

conda install -c bioconda bcftools
conda install -c bioconda bwa-mem2
conda install -c bioconda mosdepth
conda install -c bioconda gatk4