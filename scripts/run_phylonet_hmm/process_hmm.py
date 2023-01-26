#!/usr/bin/env python

import os 
import json
import operator
from operator import itemgetter
import statistics
import csv

# Directory where vcf2phylip was run
original_vcf2phylip_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/"

# Root directory where phylonet_hmm was run
root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm/"

# Tsv file with scaffold names and length in base pairs
scaffold_info_file = "scaffolds.tsv"

posterior_probability_threshold = 90

# Create dictionary containing scaffold names and their lengths in base pairs 
def create_scaffold_dict():
	with open(scaffold_info_file,"r") as f:
		scaffold_list = f.read().splitlines()

	for scaffold in scaffold_list:
		scaffold, length = scaffold.split("\t")
		scaffold_dict[scaffold] = length


# Return list of phylonet_hmm output files matched with their respective coordinate files from the scaffold alignments
def get_file_paths_pairs_list():
	
	# Create file of paths of global coordinate files for each scaffold. Save this file as "coorindate_file list" 
	os.system('find ' + original_vcf2phylip_dir + ' -name "coordinates" -type f > coordinate_file_list')
		
	# Create file of paths of rawOutput.json files for each scaffold produced by phylonet_hmm. We only want the rawOutput.json file from the "bestrun" folder. Save this file as output_file_list in hmm directory
	os.system('find ' + root_dir + ' -name "rawOutput.json" -type f | grep "best" > output_file_list')
	
	# Get zipped list matching rawOutput.json file path to global coordinate file path for each scaffold. [[Scaffold_1_rawOutput.json path, Scaffold_1_coordinates], []]
	with open("output_file_list","r") as f1:
		output_file_path_list = f1.read().splitlines()
	
	with open("coordinate_file_list","r") as f2:
		coordinate_file_path_list = f2.read().splitlines()
	
	sorted_zipped_list = [list(n) for n in zip(sorted(output_file_path_list),sorted(coordinate_file_path_list))]

	os.system('rm coordinate_file_list')
	os.system('rm output_file_list')
	
	return sorted_zipped_list

### Helper Functions for process_single_scaffold() function ###

#Get list of chromosomal coordinates from SNV alignments. Each coordinate has a corresponding posterior probability of introgression in "introgression_probabilites" list
def get_coordinate_list(coordinate_file_path):
	coordinate_list = []
	with open(coordinate_file_path,"r") as coordinate_file:
		coordinates = coordinate_file.read().splitlines()

	return coordinates

#Get list of introgression probabilities. Save to "introgression_probabilites" variable
def get_introgression_probabilities_list(json_file_path):
	with open(json_file_path,"r") as probability_file:
		introgression_probabilities = json.load(probability_file).get("posteriorProbabilityOfSpeciesTrees")[1]
	introgression_probabilities = [probability * 100 for probability in introgression_probabilities]

	return introgression_probabilities 

def get_tracts(introgression_probabilities,coordinates):
	# Create list of lists of introgression tracts in format [[start_index,stop_index,length], [start_index,stop_index,length]]
	index_tract_list = []
	coordinate_tract_list = []
	start = 0
	stop = 0
	# Find indexes of introgression tracts and add them to index_tract_list in the form of [index of start, index of stop]
	i = 0
	while (i < (len(introgression_probabilities) -1)):
		if introgression_probabilities[i] >= posterior_probability_threshold:
			start = i
			stop = i
			while introgression_probabilities[i] >= posterior_probability_threshold and i < (len(introgression_probabilities)-1):
				stop = i
				i = i+1
			index_tract_list.append([start, stop])
		i=i+1
	
	# Append each introgression tract from index_tract_list to coordinate_tract_list in format of [start_coordinate, stop_coordinate, length in base pairs] 
	for tract in index_tract_list:
		start_coordinate = coordinates[tract[0]]
		stop_coordinate = coordinates[tract[1]]
		length = int(stop_coordinate.split(":")[1]) - int(start_coordinate.split(":")[1]) + 1
		if length > 1:
			coordinate_tract_list.append([start_coordinate, stop_coordinate, length])
	
	# Sort coordinate_tract_list by index 2 of each list (length in bp) in order of highest to lowest 
	coordinate_tract_list.sort(reverse=True, key=itemgetter(2))

	return coordinate_tract_list

# Get list of all tract lengths 
def get_tract_length_dist(tract_list):
	tract_length_dist = [n[2] for n in tract_list]
	return tract_length_dist

def get_number_sites_introgressed(probability_lst):
	number_sites_introgressed = len([probability for probability in probability_lst if probability >= posterior_probability_threshold])
	return number_sites_introgressed

#### Process phylonet_hmm output for each scaffold ###

def process_single_scaffold(json_file_path,coordinate_file_path):
	coordinate_list = get_coordinate_list(coordinate_file_path)
	introgression_probabilities = get_introgression_probabilities_list(json_file_path)
	scaffold_name = str(coordinate_list[0].split(":")[0])
	
	coordinate_tract_list = get_tracts(introgression_probabilities,coordinate_list)
	
	# Get total sites on scaffold alignment
	number_sites_nexus_scaffold = len(introgression_probabilities)
	total_sites_nexus_alignments.append(number_sites_nexus_scaffold)

	# Get number of base pair sites that were declared introgressed at given threshold and append to total_snv_introgressed list. 
	number_sites_introgressed = get_number_sites_introgressed(introgression_probabilities)
	total_snv_introgressed.append(number_sites_introgressed)
	
	# Calculate percentage of sites that were declared introgressed 
	percent_sites_introgressed = ((len([probability for probability in introgression_probabilities if probability >= posterior_probability_threshold]))/len(introgression_probabilities))*100
	
	# Total length of scaffold analyzed in bases 
	total_length_scaffold_analyzed = ((int(coordinate_list[len(coordinate_list)-1].split(":")[1])) - (int(coordinate_list[0].split(":")[1])))
	total_length_nexus_alignments.append(total_length_scaffold_analyzed)

	# Actual length of scaffold
	actual_length_scaffold = int(scaffold_dict[scaffold_name])

	# Get tract length distribution 
	tract_length_dist = get_tract_length_dist(coordinate_tract_list)

	# Total number of tracts for scaffold 
	number_of_tracts = len(tract_length_dist)

	# Number of tracts on scaffold >= 10kb in length
	number_ten_kb_tracts = len([n for n in tract_length_dist if n >= 10000])

	# Append total number of tracts to total_number_tracts list 
	total_number_tracts.append(number_of_tracts)

	# Get combined length of all tracts
	combined_length_tracts = sum(tract_length_dist)

	# Calculate percentage of the scaffold that was introgressed. 
	percent_scaffold_alignment_introgressed = (combined_length_tracts/total_length_scaffold_analyzed)*100

	percent_scaffold_introgressed = (combined_length_tracts/actual_length_scaffold) * 100 

	combined_results.append(coordinate_tract_list)

	results = [scaffold_name, number_sites_nexus_scaffold, number_sites_introgressed, percent_sites_introgressed, total_length_scaffold_analyzed, actual_length_scaffold, number_of_tracts, number_ten_kb_tracts, combined_length_tracts, percent_scaffold_alignment_introgressed, percent_scaffold_introgressed]

	return results

def write_summary_stats(combined_tract_length_distribution, combined_ten_kb_tract_length_distribution):
	# Get total number of sites in concatenated scaffold alignments
	total_nexus_sites = sum(total_sites_nexus_alignments)
	print("Total nexus sites tested: {}".format(total_nexus_sites))

	# Get total length of concatenated scaffolds analyzed (The first SNV is not necesarily the first position on the scaffold)
	total_nexus_length = sum(total_length_nexus_alignments)
	print("Total nexus length tested: {}".format(total_nexus_length))

	total_actual_length = 0
	for value in scaffold_dict.values():
		total_actual_length += int(value)
	print("Total actual length of scaffolds analyzed: {}".format(total_actual_length))

	total_number_sites_introgressed = sum(total_snv_introgressed)
	print("Total number Sites Introgressed: {}".format(total_number_sites_introgressed))

	percent_snv_introgressed = (total_number_sites_introgressed/total_nexus_sites) * 100
	print("Percent SNV sites introgressed: {}".format(percent_snv_introgressed))

	total_bases_introgressed = sum(combined_tract_length_distribution)
	print("Total bases introgressed: {}".format(total_bases_introgressed))

	print("Total number of tracts: {}".format(sum(total_number_tracts)))

	total_length_introgression_tracts = sum(combined_tract_length_distribution)
	print("Total length of all introgression tracts: {}".format(total_length_introgression_tracts))

	percent_genome_analyzed_introgressed = (total_length_introgression_tracts/total_nexus_length) * 100
	print("Percent genome (analyzed) introgressed: {}".format(percent_genome_analyzed_introgressed))

	percent_genome_analyzed = (total_nexus_length / total_actual_length) * 100
	print("Percent scaffolds analyzed by Phylonet HMM: {}".format(percent_genome_analyzed))

	percent_genome_actual_introgressed = (total_length_introgression_tracts/total_actual_length) * 100
	print("Percent scaffolds (actual) introgressed: {}".format(percent_genome_actual_introgressed))

	# Calculate median length of all tracts
	median_tract_length = statistics.median(combined_tract_length_distribution)
	print("Median tract length: {}".format(median_tract_length))

	# Calculate mean length of all tracts
	mean_tract_length = statistics.mean(combined_tract_length_distribution)
	print("Mean tract length: {}".format(mean_tract_length))

	# Calculate standard deviation of all tracts
	stdev_tract_length = statistics.stdev(combined_tract_length_distribution)
	print("Standard Deviation of all tracts: {}".format(stdev_tract_length))

	num_ten_kb_tracts = len(combined_ten_kb_tract_length_distribution)
	print("Number of 10kb tracts: {}".format(num_ten_kb_tracts))

	total_length_ten_kb_introgression_tracts = sum(combined_ten_kb_tract_length_distribution)
	print("Total length of 10kb introgression tracts: {}".format(total_length_ten_kb_introgression_tracts))

	# Calculate mean, median, stdev of tracts larger than 10kb
	print("Mean tract length of tracts greater than 10 kb: {}".format(statistics.mean(combined_ten_kb_tract_length_distribution)))
	print("Median tract length of tracts longer than 10 kb: {}".format(statistics.median(combined_ten_kb_tract_length_distribution)))
	print("Standard deviation of tracts greater than 10kb: {}".format(statistics.stdev(combined_ten_kb_tract_length_distribution)))

### Helper functions to write introgression tracts and tract length distributions to bed and csv files ###

# Write tract_list in format chr start stop
def write_tracts_to_bed(file_name, tract_list):
	with open(file_name,"w") as bed_file:
		for tract in tract_list:
			scaffold = str(tract).split(":")[0].split("'")[1]
			start = str(tract).split(":")[1].split("'")[0]
			stop = str(tract[1]).split(":")[1].split("'")[0]
			bed_file.write(scaffold + "\t" + start + "\t" + stop + "\t" + scaffold + ":" + start + "_" + stop + "\n")

# Write tract legnth distribution to csv
def write_tract_dist_to_csv(file_name, tract_dist):
	with open(file_name,"a") as hist_csv:
		csv_string = ','.join([str(tract_length) for tract_length in tract_dist])
		hist_csv.write(csv_string)

def main():
	# List of coordinate_tract_lists for each scaffold in format of [[start_coordinate, stop_coordinate, length in bp], []]
	global combined_results
	combined_results = []

	# List of total number of tracts for each scaffold analyzed. 
	global total_number_tracts
	total_number_tracts = []

	# List of number of sites introgressed (prob>=90) on each scaffold 
	global total_snv_introgressed
	total_snv_introgressed = []

	# List of total number of variant sites on each scaffold alignment 
	global total_sites_nexus_alignments
	total_sites_nexus_alignments = []

	# List of the total length of each scaffold analyzed. Accounts for fact that the first site in the alignment is not necessarily the first site on the physical scaffold (same goes for last site/base)
	global total_length_nexus_alignments
	total_length_nexus_alignments = []

	# Dictionary containing scaffold names and their lengths in base pairs 
	global scaffold_dict
	scaffold_dict = {}

	# Create scaffold_dict in format of {"scaffold_name": length} using scaffold_info_file specified at the top
	create_scaffold_dict()

	# Get a list of lists of [phylonet_hmm output file, scaffold coordinate file]
	files_by_scaffold_list = get_file_paths_pairs_list()

	# Initialize a csv that will contain the phylonet_hmm results separated by scaffold analyzed
	csv_file = open("results_by_scaffold.csv","w")
	writer = csv.writer(csv_file)	
	header = ["scaffold", "SNV sites", "SNV sites introgressed", "percent sites introgressed", "scaffold length tested", "scaffold actual length", "number of introgression tracts", "number of introgression tracts >= 10kb", "combined tract length", "percent scaffold (analyzed) introgressed", "percent scaffold (actual) introgressed"]
	writer.writerow(header)

	# Process phylonet_hmm output from each scaffold, write the results to the results_by_scaffold.csv file, and append results to aggregated result arrays for further manipulation.
	for scaffold_file_pair in files_by_scaffold_list:
		json_file_path = scaffold_file_pair[0]
		coordinate_file_path = scaffold_file_pair[1]
		writer.writerow(process_single_scaffold(json_file_path,coordinate_file_path))

	csv_file.close()

	# Flatten the combined_results list of lists into a single list in format of [[start_coordinate, stop_coordinate, length in bp],[]]
	combined_results = [item for sublist in combined_results for item in sublist]
	
	# Sort flat_combined_results list of all introgression tracts by tract length in base pairs in descending order. 
	combined_results.sort(reverse=True, key=itemgetter(2))

	# Get list of introgression tracts that are larger than 10kb
	ten_kb_tracts = [n for n in combined_results if int(n[2]) >= 10000]
	
	# Get list of all tract lengths
	combined_tract_length_distribution = get_tract_length_dist(combined_results)

	# Get list of all tract lengths greater than 10kb in length
	combined_ten_kb_tract_length_distribution = [n for n in combined_tract_length_distribution if n>= 10000]

	write_summary_stats(combined_tract_length_distribution, combined_ten_kb_tract_length_distribution)

	write_tracts_to_bed("tracts.bed", combined_results)
	write_tracts_to_bed("ten_kb_tracts.bed", ten_kb_tracts)

	write_tract_dist_to_csv("tract_dist.csv", combined_ten_kb_tract_length_distribution)

if __name__ == "__main__":
        main()
