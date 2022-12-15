import os
import sys
import random
import statistics
import subprocess
import multiprocessing
from joblib import Parallel, delayed

num_cores = multiprocessing.cpu_count()

# Number of bootstrap replicates
num_replicates = 10000

root_dir = "/Users/matt/Documents/Github/dissertation/scripts/delta/"
 
bootstrap_dir = root_dir + "bootstrap/"
make_bootstrap_dir = "mkdir -p {}".format(bootstrap_dir)
os.system(make_bootstrap_dir)

tree_file = "dep-franc"
species_1 = "depressus"
species_2 = "franciscanus"
species_3 = "nudus"

def get_tree_list(file):
	return open(file,"r").read().splitlines()

def count_trees(trees):
	topology_1 = "({},({},{})".format(species_1, species_2, species_3)
	topology_2 = "({},({},{})".format(species_1, species_3, species_2)
	topology_3 = "({},({},{})".format(species_2, species_1, species_3)
	topology_4 = "({},({},{})".format(species_2, species_3, species_1)
	topology_5 = "({},({},{})".format(species_3, species_1, species_2)
	topology_6 = "({},({},{})".format(species_3, species_2, species_1)

	top_1_count = 0
	top_2_count = 0
	top_3_count = 0

	for tree in trees:
		if topology_1 in tree or topology_2 in tree:
			top_1_count += 1
		
		elif topology_3 in tree or topology_4 in tree:
			top_2_count += 1
		
		elif topology_5 in tree or topology_6 in tree:
			top_3_count += 1

	delta = (top_3_count - top_2_count) / (top_3_count + top_2_count)

	dp = (top_3_count - top_2_count) / (top_3_count + top_2_count + top_1_count)

	return delta

def get_random_tree(trees):
	random_number = random.randint(0,len(trees)-1)
	random_tree = trees[random_number]
	return random_tree

def create_pseudoreplicate(number, trees):
	with open("replicate_" + str(number),"a") as resampled_single_locus_trees_file:
		for n in range(0, len(trees)-1):
			random_tree = get_random_tree(trees)
			resampled_single_locus_trees_file.write(random_tree + "\n")

def create_delta_value_list():
	delta_value_list =[]
	replicate_file_list = os.listdir(bootstrap_dir)
	for file in replicate_file_list:
		with open(file,"r") as f:
			trees = f.read().splitlines()
			delta_value_list.append(count_trees(trees))
	return delta_value_list

def generate_output(delta_value_list, raw_delta_value):
	z_score = str((raw_delta_value)/(statistics.stdev(delta_value_list)))
	with open("bootstrap_delta_results.txt","a") as result_file:
		result_file.write("Raw Delta: " + str(raw_delta_value) + "\n")
		result_file.write("Number of bootstrap replicates: " + str(len(delta_value_list)) + "\n")
		result_file.write("Mean Delta: " + str(statistics.mean(delta_value_list)) + "\n")
		result_file.write("Standard Deviation Delta: " + str(statistics.stdev(delta_value_list)) + "\n")
		result_file.write("Z score for Delta: " + z_score)

def main():
	trees = get_tree_list(tree_file)
	raw_delta_value = count_trees(trees)
	os.chdir(bootstrap_dir)
	Parallel(n_jobs=num_cores)(delayed(create_pseudoreplicate)(i,trees) for i in range(1, num_replicates + 1))
	delta_value_list = create_delta_value_list()
	os.chdir(root_dir)
	generate_output(delta_value_list,raw_delta_value)

if __name__ == "__main__":
	main()
