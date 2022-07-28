#!/usr/bin/env python

import os 
import json
import operator
from operator import itemgetter
import statistics

original_vcf2phylip_dir = "/hb/groups/pogson_group/vcf2phylip_4way/"
working_directory = "/hb/groups/pogson_group/phylonet/4way_100runs/"

results_by_scaffold = []
combined_results = []
total_number_tracts = []
total_number_sites = []
total_sites_nexus_alignments = []
total_length_nexus_alignments = []

def get_file_paths_pairs_list():
	os.chdir(original_vcf2phylip_dir)
	os.system('find "$(pwd)" -name "coordinates" -type f > ' + working_directory + 'coordinate_file_list')
	os.chdir(working_directory)
	os.system('find "$(pwd)" -name "rawOutput.json" -type f | grep "best" > output_file_list')
	
	f1 = open("output_file_list","r")
	output_file_list = f1.read().splitlines()
	
	f2 = open("coordinate_file_list","r")
	coordinate_file_list = f2.read().splitlines()
	
	sorted_output_file_list = sorted(output_file_list)
	sorted_coordinate_file_list = sorted(coordinate_file_list)
	zipped_list = [list(n) for n in zip(sorted_output_file_list,sorted_coordinate_file_list)]

	f3 = open("coordinates_and_json_file","w")
	for n in zipped_list:
		f3.write(str(n)+"\n")
	f3.close()

	file_paths_pairs = open("coordinates_and_json_file","r")
	file_paths_pairs_lines = file_paths_pairs.read().splitlines()

	file_paths_pairs_list = []

	for n in file_paths_pairs_lines:
		output_file = n.split(",")[0][1:].split("'")[1]
		coordinate_file = n.split(",")[1].split("]")[0].split("'")[1]
		file_paths_pairs_list.append([output_file,coordinate_file])

	return file_paths_pairs_list

def find_tracts(introgression_probabilities,coordinates):
	index_tract_list = []
	start = 0
	stop = 0
	# Find indexes of introgression tracts and add them to index_tract_list in the form of [index of start, index of stop]
	i = 0
	while (i < (len(introgression_probabilities) -1)):
		if introgression_probabilities[i] >= 90:
			start = i
			stop = i
			while introgression_probabilities[i] >= 90 and i < (len(introgression_probabilities)-1):
				stop = i
				i = i+1
			index_tract_list.append([start, stop])
		i=i+1
	# Add length in base pairs to index_tract_list (indexes of introgression tracts)
	for n in index_tract_list:
		x = n[0]
		y = n[1]
		x2 = coordinates[x]
		y2 = coordinates[y]
		z = int(y2.split(":")[1]) - int(x2.split(":")[1]) + 1
		n.append(z)
	return index_tract_list

def get_coordinate_tract_list_and_tract_lengths(index_tract_list,coordinates):
	# List of the [start,stop] of all introgression tracts from tract_list
	coordinate_tract_list = []
	# Open coordinate file. Each line corresponds to the chromosomal coordinate of the matching index in the introgression_probabilities list
	for n in index_tract_list:
		x = n[0]
		y = n[1]
		x2 = coordinates[x]
		y2 = coordinates[y]
		coordinate_tract_list.append([x2,y2])
	# Add length in base pairs to each introgression tract 
	for n in coordinate_tract_list:
		z = int(n[1].split(":")[1]) - int(n[0].split(":")[1]) + 1
		n.append(z)
	return coordinate_tract_list

def sort_index_tract_list(index_tract_list):
	sorted_index_tract_list = (sorted(index_tract_list, key=itemgetter(2), reverse=True))
	return sorted_index_tract_list

def sort_coordinate_tract_list(coordinate_tract_list):
	sorted_coordinate_tract_list = (sorted(coordinate_tract_list, key=itemgetter(2), reverse=True))
	return sorted_coordinate_tract_list

def get_tract_length_dist(sorted_coordinate_tract_list):
	tract_length_dist = [n[2] for n in sorted_coordinate_tract_list]
	return tract_length_dist

def process_single_scaffold(json_file_path,coordinate_file_path):
	json_file = open(json_file_path,"r")
	coordinate_file = open(coordinate_file_path,"r")

	# List of chromosomal coordinates. Each coordinate has a corresponding posterior probability of introgression in "introgression_probabilites" list
	coordinates = []
	inputs = coordinate_file.read().splitlines()
	# Reformat coordinates into a list
	[coordinates.append(item)for item in inputs]

	scaffold_name = []
	scaffold_name.append(coordinates[0].split(":")[0])
	scaffold_name = str(scaffold_name[0])

	# Get list of introgression probabilities. Save to "introgression_probabilites" variable
	probability_file =open(json_file_path,'r')
	introgression_probabilities=json.load(probability_file).get("posteriorProbabilityOfSpeciesTrees")[1]
	introgression_probabilities = [probability * 100 for probability in introgression_probabilities]
	length = len(introgression_probabilities)

	index_tract_list = find_tracts(introgression_probabilities,coordinates)
	coordinate_tract_list = get_coordinate_tract_list_and_tract_lengths(index_tract_list,coordinates)
	sorted_index_tract_list = sort_index_tract_list(index_tract_list)
	sorted_coordinate_tract_list = sort_coordinate_tract_list(coordinate_tract_list)
	tract_length_dist = get_tract_length_dist(sorted_coordinate_tract_list)

	number_sites_introgressed = len([probability for probability in introgression_probabilities if probability >= 90])
	#print(number_sites_introgressed)
	total_number_sites.append(number_sites_introgressed)
	
	sum_tract_length_dist = sum(tract_length_dist)
	percent_sites_introgressed = ((len([probability for probability in introgression_probabilities if probability >= 90]))/len(introgression_probabilities))*100
	percent_scaffold_alignment_introgressed = ((sum(tract_length_dist))/((int(coordinates[len(coordinates)-1].split(":")[1])) - (int(coordinates[0].split(":")[1]))))*100
	
	# Get total sites on scaffold alignment
	number_sites_nexus_scaffold = len(introgression_probabilities)
	total_sites_nexus_alignments.append(number_sites_nexus_scaffold)
	
	# Get total length in bases of nexus scaffold alignment based of genomic coordinates
	total_length_scaffold_analyzed = ((int(coordinates[len(coordinates)-1].split(":")[1])) - (int(coordinates[0].split(":")[1])))
	total_length_nexus_alignments.append(total_length_scaffold_analyzed)

	#print(scaffold_name)
	#print(percent_sites_introgressed)
	#print(percent_scaffold_alignment_introgressed)
	#print(scaffold_name,percent_scaffold_alignment_introgressed,percent_sites_introgressed,percent_scaffold_alignment_introgressed - percent_sites_introgressed)

	results_by_scaffold.append([scaffold_name,sorted_coordinate_tract_list,sorted_index_tract_list,tract_length_dist])

	combined_results.append(coordinate_tract_list)

	number_of_tracts = len(tract_length_dist)
	total_number_tracts.append(number_of_tracts)
	
def main():
	
	for pair in get_file_paths_pairs_list():
		json_file_path = pair[0]
		coordinate_file_path = pair[1]
		process_single_scaffold(json_file_path,coordinate_file_path)

	flat_combined_results = [item for sublist in combined_results for item in sublist]
	sorted_combined_results = sort_coordinate_tract_list(flat_combined_results)
	
	# Get total number of sites in all nexus scaffolds analzyed
	#print("Total nexus sites tested:")
	total_nexus_sites = sum(total_sites_nexus_alignments)
	#print(total_nexus_sites)
	
	# Get total length of scaffolds analyzed
	#print("Total nexus length tested:")
	total_nexus_length = sum(total_length_nexus_alignments)
	#print(total_nexus_length)

	#print(sorted_combined_results[0:10])

	print("Total Number Sites Introgressed:")
	print(sum(total_number_sites))
	num_ten_kb_tracts = len([n for n in sorted_combined_results if int(n[2]) >= 10000])
	#print("Number of ten kb tracts:")
	#print(num_ten_kb_tracts)
	
	#print("Total Number of tracts:")
	#print(len(sorted_combined_results))
	#print(sorted_combined_results[0:10])
	combined_tract_length_distribution = get_tract_length_dist(sorted_combined_results)
	
	# Print combined tract length distribution
	
	print("Total bases introgressed")
	print(sum(combined_tract_length_distribution))	


	combined_ten_kb_or_greater_tract_length_distribution = [n for n in combined_tract_length_distribution if n>= 10000]
	#print("Number of ten kb tracts:")
	#print(len(combined_ten_kb_or_greater_tract_length_distribution))
	sorted_combined_ten_kb_or_greater_tract_length_distribution = sorted(combined_ten_kb_or_greater_tract_length_distribution,reverse=True)
	
	# Calculate median length of tracts larger than 10kb
	#print("Mean tract length for tracts greater than 10 kb:")
	#print(statistics.mean(sorted_combined_ten_kb_or_greater_tract_length_distribution))
	#print("Median tract length for tracts longer than 10 kb:")
	#print(statistics.median(sorted_combined_ten_kb_or_greater_tract_length_distribution))
	#print("Standard deviation for tracts greater than 10kb:")
	#print(statistics.stdev(sorted_combined_ten_kb_or_greater_tract_length_distribution))
	
	# Calculate median length of all tracts
	median_tract_length = statistics.median(combined_tract_length_distribution)
	#print("median length of all tracts:")
	#print(median_tract_length)
	#mean_tract_length = statistics.mean(combined_tract_length_distribution)
	#print("Mean length of all tracts:")
	#print(mean_tract_length)
	#print("Standard Deviation of all tracts:")
	#print(statistics.stdev(combined_tract_length_distribution))
	#print(combined_tract_length_distribution)
	
	# Total number of tracts
	#print("Total number of tracts:")
	#print(sum(total_number_tracts))

	# Need to write these results to files now that can be easily parsed.
	
	# Write sorted_combined_results to combined_coordinate_tract_list.bed in format chr start stop
	bed_file = open("tracts.bed","w")
	for n in sorted_combined_results:
		chr = str(n).split(":")[0].split("'")[1]
		start = str(n).split(":")[1].split("'")[0]
		stop = str(n[1]).split(":")[1].split("'")[0]
		bed_file.write(chr + "\t" + start + "\t" + stop + "\t" + chr + ":" + start + "_" + stop + "\n")
	bed_file.close()

	# Write combined_tract_length_distribution to csv
	hist_csv = open("tract_dist.csv","a")
	for n in sorted_combined_ten_kb_or_greater_tract_length_distribution:
		hist_csv.write(str(n) + ",")
	hist_csv.close()

	# Write results_by_scaffold to json file object???

if __name__ == "__main__":
        main()
