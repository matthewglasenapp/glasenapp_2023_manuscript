Before running run_phylonet_hmm.sh, must update each python script with relevant info. 

1) scaffolds.py creates nexus alignments for each scaffold as well as a coordinate file reporting the global coordinate of each SNV site applied in the nexus alignment. It first subsets a multisample vcf by scaffold and stores the output in root_dir/vcf_by_scaffold/. It then creates the nexus matrixes using the vcf2phylip program and stores the output in root_dir/scaffold_nexus_alignments/[scaffold]/.

2) create_matrices.py creates the phylonet_hmm nexus input file for each scaffold using the output from the vcf2phylip program. The phylonet_hmm input files are store at /hb/scratch/mglasena/phylonet_hmm_input/hmm_nexus_files/

3) run_phylonet_hmm.py runs Phylonet HMM on each scaffold alignment input nexus file. 

4) process_hmm.py processes the output files of Phylonet HMM into useful stats. scaffolds.tsv - a tsv file with scaffold names and length in base pairs, must be in the directory where run_phylonet_hmm.sh is launched. 