import os
import gzip

# /hb/scratch/mglasena/run_mosdepth/nudus_SRR5767281.regions.bed.gz

bed_file = "/hb/scratch/mglasena/run_mosdepth/protein_coding_genes.bed"

bed_file_dir = "/hb/scratch/mglasena/run_mosdepth/"

min_cov_threshold = 10
max_cov_threshold = 150

gene_dict = dict()

new_dict = dict()

def get_bed_file_list():
	#get_bed_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > region_files".format(bed_file_dir)
	#os.system(get_bed_file_paths)

	with open("region_files", "r") as f:
		file_list = f.read().splitlines()

	return file_list

def create_gene_dict():
	with open(bed_file,"r") as f:
		for line in f:
			gene = line.split("\t")[3]
			gene_dict[gene] = []

def filter_genes(file):
	with gzip.open(file,"rt") as f:
		for line in f:
			gene_id = line.split("\t")[-2]
			mean_depth = float(line.split("\t")[-1].strip())
			gene_dict[gene_id].append(mean_depth) 

def main():
	bed_file_list = get_bed_file_list()
	
	create_gene_dict()

	for file in bed_file_list:
		filter_genes(file)

	for key, value in gene_dict.items():
		if min(value) > min_cov_threshold and max(value) < max_cov_threshold:
			new_dict[key] = value

	print(len(new_dict))
	for key,value in new_dict.items():
		print("{} : {}".format(key,value))
	
	#print(new_dict)

if __name__ == "__main__":
	main()


