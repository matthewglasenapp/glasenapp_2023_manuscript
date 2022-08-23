import os
from joblib import Parallel, delayed

bam_file_dir = "/hb/scratch/mglasena/data/dedup_mapped_bam_files/"

def get_file_paths_list():
	get_bam_paths_file = 'find {} -type f -name *.bam* | grep -v "bai" > bam_files.txt'.format(bam_file_dir)
	os.system(get_bam_paths_file)

	with open("bam_files.txt", "r") as f:
		bam_file_paths_list = f.read().splitlines()

	return bam_file_paths_list

def run_flagstat(bam_file):
	output_file = bam_file.split(".")[0] + "_primary_flagstat.tsv"
	flagstat = "samtools view -u -F 256 {} | samtools flagstat -O tsv - > {}".format(bam_file, output_file)
	os.system(flagstat)

def main():
	bam_file_paths_list = get_file_paths_list()

	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_flagstat)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()