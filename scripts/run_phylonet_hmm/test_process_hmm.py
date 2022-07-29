#!/usr/bin/env python

import os 
import json
import operator
from operator import itemgetter
import statistics
#import multiprocessing
#from joblib import Parallel, delayed
#num_cores = multiprocessing.cpu_count()

#Directory where vcf2phylip was run
original_vcf2phylip_dir = "/hb/groups/pogson_group/vcf2phylip/vcf2phylip_4way/"

#Root directory where phylonet_hmm was run
root_dir = "/hb/groups/pogson_group/phylonet/4way_100runs/"

posterior_probability_threshold = 90

# List of lists of results for each scaffold in format of [[scaffold_name,sorted_coordinate_tract_list,sorted_index_tract_list,tract_length_dist]]
results_by_scaffold = []

# List of coordinate_tract_lists for each scaffold in format of [[start_coordinate, stop_coordinate, length in bp], []]
combined_results = []

# List of total number of tracts for each scaffold analyzed. 
total_number_tracts = []

# List of number of sites introgressed on each scaffold 
total_number_sites = []

# List of total number of variant sites on each scaffold alignment 
total_sites_nexus_alignments = []

# List of the total length of each scaffold analyzed. Accounts for fact that the first site in the alignment is not necessarily the first site on the physical scaffold (same goes for last site/base)
total_length_nexus_alignments = []

# Return list of 
def get_file_paths_pairs_list():
	
	#Create file of paths of global coordinate files for each scaffold. Save this file as "coorindate_file list" 
	os.system('find ' + original_vcf2phylip_dir + ' "$(pwd)" -name "coordinates" -type f > coordinate_file_list')
		
	#Create file of paths of rawOutput.json files for each scaffold produced by phylonet_hmm. We only want the rawOutput.json file from the "bestrun" folder. Save this file as output_file_list in hmm directory
	os.system('find ' + root_dir + ' "$(pwd)" -name "rawOutput.json" -type f | grep "best" > output_file_list')
	
	#Get zipped list matching rawOutput.json file path to global coordinate file path for each scaffold. [[Scaffold_1_rawOutput.json path, Scaffold_1_coordinates], []]
	with open("output_file_list","r") as f1:
		output_file_path_list = f1.read().splitlines()
	
	with open("coordinate_file_list","r") as f2:
		coordinate_file_path_list = f2.read().splitlines()
	
	sorted_zipped_list = [list(n) for n in zip(sorted(output_file_path_list),sorted(coordinate_file_path_list))]

	os.system('rm coordinate_file_list')
	os.system('rm output_file_list')
	
	return sorted_zipped_list

#Get list of chromosomal coordinates from SNV alignments. Each coordinate has a corresponding posterior probability of introgression in "introgression_probabilites" list
def get_coordinate_list(coordinate_file_path):
	coordinate_list = []
	with open(coordinate_file_path,"r") as coordinate_file:
		coordinates = coordinate_file.read().splitlines()
	[coordinate_list.append(coordinate)for coordinate in coordinates]

	return coordinate_list

#Get list of introgression probabilities. Save to "introgression_probabilites" variable
def get_introgression_probabilities_list(json_file_path):
	with open(json_file_path,"r") as probability_file:
		introgression_probabilities = json.load(probability_file).get("posteriorProbabilityOfSpeciesTrees")[1]
	introgression_probabilities = [probability * 100 for probability in introgression_probabilities]

	return introgression_probabilities 

def find_tracts(introgression_probabilities,coordinates):
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
	
	# Append length in base pairs to each tract in index_tract_list (indexes of introgression tracts)
	for tract in index_tract_list:
		start_coordinate = coordinates[tract[0]]
		stop_coordinate = coordinates[tract[1]]
		length = int(stop_coordinate.split(":")[1]) - int(start_coordinate.split(":")[1]) + 1
		tract.append(length)
		coordinate_tract_list.append([start_coordinate, stop_coordinate, length])
	
	return index_tract_list, coordinate_tract_list

# Sort index_tract_list by index 2 of each list (length in bp) in order of highest to lowest 
def sort_tract_list(tract_list):
	sorted_tract_list = (sorted(tract_list, key=itemgetter(2), reverse=True))
	return sorted_tract_list

# Get list of all tract lengths 
def get_tract_length_dist(tract_list):
	tract_length_dist = [n[2] for n in tract_list]
	return tract_length_dist

def process_single_scaffold(json_file_path,coordinate_file_path):
	coordinate_list = get_coordinate_list(coordinate_file_path)

	# Get name of scaffold
	scaffold_name = str(coordinate_list[0].split(":")[0])
	
	introgression_probabilities = get_introgression_probabilities_list(json_file_path)

	find_tracts_results = find_tracts(introgression_probabilities,coordinate_list)
	index_tract_list = find_tracts_results[0]
	coordinate_tract_list = find_tracts_results[1]
	
	# Sort index_tract_list and coordinate_tract_list by index 2 of each list (length in bp) in order of highest to lowest 
	sorted_index_tract_list = sort_tract_list(index_tract_list)
	sorted_coordinate_tract_list = sort_tract_list(coordinate_tract_list)
	
	# Get list of all tract lengths 
	tract_length_dist = get_tract_length_dist(sorted_coordinate_tract_list)

	# Get number of base pair sites that were declared introgressed at given threshold and append to total_number_sites list. 
	number_sites_introgressed = len([probability for probability in introgression_probabilities if probability >= posterior_probability_threshold])
	#print(number_sites_introgressed)
	total_number_sites.append(number_sites_introgressed)
	
	# Get sum of tract_length_dist list, which is a list of all tract lengths on scaffold
	sum_tract_length_dist = sum(tract_length_dist)
	
	# Calculate percentage of sites that were declared introgressed 
	percent_sites_introgressed = ((len([probability for probability in introgression_probabilities if probability >= posterior_probability_threshold]))/len(introgression_probabilities))*100
	
	# Calculate percentage of the scaffold that was introgressed. 
	percent_scaffold_alignment_introgressed = ((sum(tract_length_dist))/((int(coordinate_list[len(coordinate_list)-1].split(":")[1])) - (int(coordinate_list[0].split(":")[1]))))*100
	
	# Get total sites on scaffold alignment
	number_sites_nexus_scaffold = len(introgression_probabilities)
	total_sites_nexus_alignments.append(number_sites_nexus_scaffold)
	
	# Get total length in bases of nexus scaffold alignment based of genomic coordinates
	total_length_scaffold_analyzed = ((int(coordinate_list[len(coordinate_list)-1].split(":")[1])) - (int(coordinate_list[0].split(":")[1])))
	total_length_nexus_alignments.append(total_length_scaffold_analyzed)

	#print(scaffold_name)
	#print(percent_sites_introgressed)
	#print(percent_scaffold_alignment_introgressed)
	#print(scaffold_name,percent_scaffold_alignment_introgressed,percent_sites_introgressed,percent_scaffold_alignment_introgressed - percent_sites_introgressed)

	results_by_scaffold.append([scaffold_name,sorted_coordinate_tract_list,sorted_index_tract_list,tract_length_dist])

	combined_results.append(coordinate_tract_list)

	# Calculate total number of tracts for scaffold 
	number_of_tracts = len(tract_length_dist)
	
	# Append total number of tracts to total_number_tracts list 
	total_number_tracts.append(number_of_tracts)

# Write sorted_flattened_combined_results to combined_coordinate_tract_list.bed in format chr start stop
def write_tracts_to_bed(tract_list):
	with open("tracts.bed","w") as bed_file:
		for tract in tract_list:
			scaffold = str(tract).split(":")[0].split("'")[1]
			start = str(tract).split(":")[1].split("'")[0]
			stop = str(tract[1]).split(":")[1].split("'")[0]
			bed_file.write(scaffold + "\t" + start + "\t" + stop + "\t" + scaffold + ":" + start + "_" + stop + "\n")

# Write combined_tract_length_distribution to csv
def write_tract_dist_to_csv(tract_dist):
	with open("tract_dist.csv","a") as hist_csv:
		for tract in tract_dist:
			hist_csv.write(str(tract) + ",")

def main():
	#Get a list of lists of [output_file,coordinate_file]
	files_by_scaffold_list = get_file_paths_pairs_list()

	# Process each scaffold and append results to aggregated result arrays
	for scaffold_file_pair in files_by_scaffold_list:
		json_file_path = scaffold_file_pair[0]
		coordinate_file_path = scaffold_file_pair[1]
		process_single_scaffold(json_file_path,coordinate_file_path)

	# Flatten the combined_results list of lists into a single list in format of [[start_coordinate, stop_coordinate, length in bp],[]]
	flat_combined_results = [item for sublist in combined_results for item in sublist]
	
	# Sort flat_combined_results list of all introgression tracts by tract length in base pairs in descending order. 
	sorted_flattened_combined_results = sort_tract_list(flat_combined_results)
	
	# Get total number of sites in all nexus scaffolds analzyed
	total_nexus_sites = sum(total_sites_nexus_alignments)
	#print("Total nexus sites tested: {}".format(total_nexus_sites))
	
	# Get total length of scaffolds analyzed
	total_nexus_length = sum(total_length_nexus_alignments)
	#print("Total nexus length tested: {}".format(total_nexus_length))

	#print(sorted_flattened_combined_results[0:10])

	print("Total Number Sites Introgressed: {}".format(sum(total_number_sites)))
	
	num_ten_kb_tracts = len([n for n in sorted_flattened_combined_results if int(n[2]) >= 10000])
	#print("Number of ten kb tracts: {}".format(num_ten_kb_tracts))
	
	#print("Total Number of tracts: {}".format(len(sorted_flattened_combined_results)))

	#print(sorted_flattened_combined_results[0:10])
	
	# Get list of all tract lengths
	combined_tract_length_distribution = get_tract_length_dist(sorted_flattened_combined_results)
	
	# Print combined tract length distribution
	print("Total bases introgressed: {}".format(sum(combined_tract_length_distribution)))

	combined_ten_kb_or_greater_tract_length_distribution = [n for n in combined_tract_length_distribution if n>= 10000]
	
	#print("Number of ten kb tracts: {}".format(len(combined_ten_kb_or_greater_tract_length_distribution)))
	
	# Get list of all tracts > 10kb in length sorted in descending order of tract length 
	sorted_combined_ten_kb_or_greater_tract_length_distribution = sorted(combined_ten_kb_or_greater_tract_length_distribution,reverse=True)
	
	# Calculate median length of tracts larger than 10kb
	#print("Mean tract length for tracts greater than 10 kb: {}".format(statistics.mean(sorted_combined_ten_kb_or_greater_tract_length_distribution)))
	
	#print("Median tract length for tracts longer than 10 kb: {}".format(statistics.median(sorted_combined_ten_kb_or_greater_tract_length_distribution)))
	
	#print("Standard deviation for tracts greater than 10kb: {}".format(statistics.stdev(sorted_combined_ten_kb_or_greater_tract_length_distribution)))
	
	# Calculate median length of all tracts
	median_tract_length = statistics.median(combined_tract_length_distribution)
	#print("median length of all tracts: {}".format(median_tract_length))

	# Calculate mean length of all tracts
	mean_tract_length = statistics.mean(combined_tract_length_distribution)
	#print("Mean length of all tracts: {}".format(mean_tract_length))

	# Calculate standard deviation of all tracts
	#print("Standard Deviation of all tracts: {}".format(statistics.stdev(combined_tract_length_distribution)))

	#print(combined_tract_length_distribution)
	
	# Total number of tracts
	#print("Total number of tracts: {}".format(sum(totoal_number_tracts)))

	write_tracts_to_bed(sorted_flattened_combined_results)

	write_tract_dist_to_csv(sorted_combined_ten_kb_or_greater_tract_length_distribution)

	# Write results_by_scaffold to json file object???

if __name__ == "__main__":
        main()
