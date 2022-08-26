# Create conda environment for analysis

module load miniconda3.9

conda create -y -n analysis

conda activate analysis

conda install -c bioconda bcftools
conda install -c bioconda bwa-mem2
conda install -c bioconda mosdepth
conda install -c bioconda gatk4


# Common conda commands
# conda info --envs
# conda remove --name env --all

# QuiBL Install
#module load miniconda3.9
#conda create -y -n quibl python=2.7
#conda activate quibl
#pip install -r quibl_requirements.txt --user
