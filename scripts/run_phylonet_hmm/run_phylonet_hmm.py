import os
import multiprocessing
from joblib import Parallel, delayed

num_jobs = 1
memory = 125

phylonet_hmm_path = "/hb/groups/pogson_group/dissertation/software/phylonet_hmm/PHiMM.jar"

input_file_dir = "/hb/scratch/mglasena/test_phylonet_hmm/hmm_input/hmm_nexus_files"

root_dir = "/hb/scratch/mglasena/test_phylonet_hmm/"
hmm_dir = root_dir + "hmm/"
os.mkdir(hmm_dir)

def get_scaffold_input_nexus_file_path_list():
	create_scaffold_input_nexus_file_paths_file = 'find {} "$(pwd)" -name "*.nexus" -type f > scaffold_input_nexus_file_paths_file'.format(input_file_dir)
	os.system(create_scaffold_input_nexus_file_paths_file)

	with open("scaffold_input_nexus_file_paths_file", "r") as f:
		scaffold_input_nexus_file_path_list = f.read().splitlines()
	
	return scaffold_input_nexus_file_path_list

def run_hmm(scaffold):
	run_hmm = "java -Xmx{}g -jar {} {}".format(memory, phylonet_hmm_path, scaffold)
	os.system(run_hmm)

def main():
	os.chdir(hmm_dir)

	scaffold_input_nexus_file_path_list = get_scaffold_input_nexus_file_path_list()
	
	#Parallel(n_jobs=num_jobs)(delayed(run_hmm)(scaffold_input_nexus_file) for scaffold_input_nexus_file in scaffold_input_nexus_file_path_list)
	for file in scaffold_input_nexus_file_path_list:
		run_hmm(file)

if __name__ == "__main__":
        main()
