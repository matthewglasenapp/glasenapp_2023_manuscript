# Phylonet-HMM Analysis Workflow

Before running run_phylonet_hmm.sh, you must update each python script with relevant info. 

1) scaffolds.py 
Create nexus alignments for each scaffold as well as a coordinate file reporting the global coordinate of each SNV site applied in the nexus alignment. It first subsets a multisample vcf by scaffold and stores the output in root_dir/vcf_by_scaffold/. It then creates the nexus matrices using the vcf2phylip program and stores the output in root_dir/scaffold_nexus_alignments/[scaffold]/.

2) create_matrices.py 
Create the phylonet_hmm nexus input file for each scaffold using the output from the vcf2phylip program. The phylonet_hmm input files are stored at /hb/scratch/mglasena/phylonet_hmm_input/hmm_nexus_files/

3) run_phylonet_hmm.py 
Run Phylonet HMM on each scaffold alignment input nexus file. 

4) process_hmm.py 
Processes the output files of Phylonet HMM into useful stats. scaffolds.tsv - a tsv file with scaffold names and length in base pairs, must be in the directory where run_phylonet_hmm.sh is launched. 

5) investigate_tracts.py

# process_hmm.py 

At the beginning of the script, you must specify the following variables:

1) original_vcf2phylip_dir: Directory where vcf2phylip was run. This script assumes that in the original_vcf2phylip_dir there is a directory for each scaffold that contains a file called "coordinates" with the applied coordinates in the scaffold alignment in the format of scaffold:pos.

2) root_dir: Root directory where phylonet_hmm was run

3) scaffold_info_file: Tsv file with scaffold names and length in base pairs

4) posterior_probability_threshold: Threshold cutoff to use for declaring sites introgressed (default=90)

## Walkthrough of main() in process_hmm.py

1) Instantiate several global list variables. These list will be continuously updated as each scaffold is processed and will be used to calculate aggregate statistics for all scaffolds at the end. 

2) Call create_scaffold_dict() to create a dictionary of scaffold names and their lengths in the format of {"scaffold_name": length}. This dictionary is created from a tsv file specified at the top of the script in the scaffold_info_file variable. The dictionary produced is used for calculations in the process_single_scaffold() function called later. 

3) Create files_by_scaffold_list variable by calling get_file_paths_pairs_list(), which creates of list of lists of [phylonet_hmm output file for scaffold, coordinate file for sites used in scaffold alignment]. This function requires specifiying the directory where vcf2phylip was run when producing scaffold SNV alignment files (original_vcf2phylip_dir) and the directory where the phylonet_hmm program was run (root_dir). Each list of file pairs will be passed to process_single_scaffold() for processing. The phylonet_hmm output file contains an array with probability of introgression for each SNV site in the scaffold alignment. The coordinate file contains the global coordinates for each site in the scaffold SNV alignment. get_file_paths_pairs_list() assumes that in the original_vcf2phylip_dir there is a directory for each scaffold that contains a file called "coordinates" with the applied coordinates in the scaffold alignment in the format of scaffold:pos. The function will extract the "best_run" file for each scaffold from the directory where phylonet_hmm was run. 

4) A csv file called "results_by_scaffold.csv" is created. This file will store the results of processing the phylonet_hmm run by scaffold. The csv file has the following header ["scaffold", "SNV sites", "SNV sites introgressed", "percent sites introgressed", "scaffold length tested", "scaffold actual length", "number of introgression tracts", "number of introgression tracts >= 10kb", combined tract length", "percent scaffold (analyzed) introgressed", "percent scaffold (actual) introgressed"]. 

5) process_single_scaffold() is called to process the phylonet_hmm output for each scaffold. The function returns a list of variables calculated while processing the phylonet_hmm output - [scaffold_name, number_sites_nexus_scaffold, number_sites_introgressed, percent_sites_introgressed, total_length_scaffold_analyzed, actual_length_scaffold, number_of_tracts, number_ten_kb_tracts, combined_length_tracts, percent_scaffold_alignment_introgressed, percent_scaffold_introgressed]. This list is written to the "results_by_scaffold.csv" file.

      process_single_scaffold() calls several other helper functions:
  
      a) get_coordinate_list() returns a list of all coordinates applied to scaffold SNV alignment 
    
      b) get_introgression_probabilities_list() returns a list the probability each site in the scaffold SNV alignment is introgressed
    
      c) get_tracts() finds introgression tracts, stretches of consecutive SNV sites with posterior probabilities of introgression exceeding the posterior_probability_threshold at the beginning of the script (default=90). This function returns a list of tracts > 1 SNV site in length
    
      d) get_number_sites_introgressed() calculates the number of SNV sites that exceeded the probability threshold for being declared introgressed
    
      e) get_tract_length_dist() returns an array of the tract lengths contained in the coordinate_tract_list returned by get_tracts()

    process_single_scaffold() also calculates numerous values and appends them to aggregate lists containing data for each scaffold:
      
      a) Calculates and appends the total number of SNV sites included in the scaffold to the total_sites_nexus_alignments list
    
      b) Calculates and appends the number of SNV sites that were declared introgressed the given probability threshold to the total_snv_introgressed list.
    
      c) Calculates the percent of SNV sites on the scaffold alignment that were declared introgressed
    
      d) Calculates and appends the total length of the analyzed scaffold alignment to the total_length_nexus_alignments list.
      
      e) Reports the actual length of the scaffold analyzed based on the current genome assembly
    
      f) Calculates and appends the total number of introgression tracts to the total_number_tracts list.
      
      g) Calculates the number of introgression tracts >= 10kb in length
      
      h) Calculates the combined length of all tracts, the percent of the analyzed scaffold declared introgressed, 
      
      i) Appends the sorted coordinate_tract_list returned by get_tracts() to the combined_results list

6) Flatten and sort the list of lists of introgression tracts by scaffold

7) Get list of introgression tracts that are larger than 10kb (ten_kb_tracts)

8) Call get_tract_length_distribution() on combined_results to get list of all tract lengths. Save to combined_tract_length_distribution

9) Save list of all tract lengths greater than 10kb in length to combined_ten_kb_tract_length_distribution

10) Use write_summary_stats() function to calculate and print out some summary stats for the aggregated data. 

11) Use write_tracts_to_bed() to create two bed files containing (i) all introgression tracts, and (ii) introgression tracts >= 10kb in length in the format of chr start stop tract_name

12) Use write_tract_dist_to_csv() to create a csv file containing a list of all tract lengths >= 10kb
