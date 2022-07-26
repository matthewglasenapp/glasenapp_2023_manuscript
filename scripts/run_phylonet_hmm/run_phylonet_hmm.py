import os
import multiprocessing
from joblib import Parallel, delayed

num_jobs = 21
memory = 5

scaffold_input_nexus_file_paths = "nexus_list.txt"

def get_scaffold_input_nexus_file_path_list():
	with open(scaffold_input_nexus_file_paths, "r") as f:
		scaffold_input_nexus_file_path_list = f.read().splitlines()
	return scaffold_input_nexus_file_path_list


def run_hmm(scaffold):
	run_hmm = "java -Xmx{}g -jar PHiMM.jar {}".format(memory, scaffold)
	os.system(run_hmm)

def main():
	scaffold_input_nexus_file_path_list = get_scaffold_input_nexus_file_path_list()
	Parallel(n_jobs=num_jobs)(delayed(run_hmm)(scaffold_input_nexus_file) for scaffold_input_nexus_file in scaffold_input_nexus_file_path_list)

if __name__ == "__main__":
        main()
