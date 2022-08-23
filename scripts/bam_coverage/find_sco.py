import os
import gzip
from itertools import islice
import csv

bed_file = "/hb/scratch/mglasena/run_mosdepth/protein_coding_genes.bed"

bed_file_dir = "/hb/scratch/mglasena/run_mosdepth/"

min_cov_threshold = 10
max_cov_threshold = 100
purpuratus_7211988_max = 150
prop_1x_threshold = 0.5

gene_dict = dict()

passed_genes_dict = dict()

def get_zipped_bed_file_list():
	#get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > regions_files".format(bed_file_dir)
	#os.system(get_regions_file_paths)

	#get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* | grep -v 'csi' > thresholds_files".format(bed_file_dir)
	#os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	return list(file_list)

def initialize_gene_dict():
	with open(bed_file,"r") as f:
		for line in f:
			gene = line.split("\t")[3]
			gene_dict[gene] = [],[]

def fill_gene_dict(regions_file, thresholds_file):
	with gzip.open(regions_file, "rt") as file1, gzip.open(thresholds_file, "rt") as file2:
		for line_file1, line_file2 in zip(file1, islice(file2, 1, None)):
			gene_id = line_file1.split("\t")[-2]
			mean_depth = float(line_file1.split("\t")[-1].strip())
			
			gene_dict[gene_id][0].append(mean_depth)

			total_base_count = int(line_file2.split("\t")[2]) - int(line_file2.split("\t")[1]) + 1
			one_x_count = int(line_file2.split("\t")[4])
			try:
				prop_1x = float(one_x_count / total_base_count)
			except ZeroDivisionError:
				prop_1x = 0.0

			gene_dict[gene_id][1].append(prop_1x)

def filter_gene_dict():
	for key, value in gene_dict.items():
		if min(value[0]) >= min_cov_threshold and max(value[0]) < max_cov_threshold and min(value[1]) >= prop_1x_threshold:
			passed_genes_dict[key] = value

	print("{} genes passed filter".format(len(passed_genes_dict)))


def write_genes_passed_filter_bed():
	with open("genes_pass_filter.bed", "a") as f1, open(bed_file,"r") as f2:
		genes_pass_filter_list = list(passed_genes_dict.keys())
		all_gene_dict = {line.split("\t")[3]:line for line in f2.read().splitlines()}

		for gene in genes_pass_filter_list:
			if gene in all_gene_dict.keys():
				f1.write(all_gene_dict[gene] + "\n")

def write_all_gene_dict_csv():
	csv_file = open("all_genes.csv","w")
	writer = csv.writer(csv_file)	
	header = ["gene", "depressus", "droebachiensis", "fragilis", "franciscanus", "intermedius", "nudus", "pallidus", "pulcherrimus_DRR107784", "pulcherrimus_SRR5767283", "purpuratus_SRR6281818", "purpuratus_SRR7211988", "depressus", "droebachiensis", "fragilis", "franciscanus", "intermedius", "nudus", "pallidus", "pulcherrimus_DRR107784", "pulcherrimus_SRR5767283", "purpuratus_SRR6281818", "purpuratus_SRR7211988"]
	writer.writerow(header)

	for key,value in gene_dict.items():
		row = [key] + value[0] + value[1]
		writer.writerow(row)

	csv_file.close()

def write_passed_genes_dict_csv():
	csv_file = open("passed_genes.csv","w")
	writer = csv.writer(csv_file)	
	header = ["gene", "depressus", "droebachiensis", "fragilis", "franciscanus", "intermedius", "nudus", "pallidus", "pulcherrimus_DRR107784", "pulcherrimus_SRR5767283", "purpuratus_SRR6281818", "purpuratus_SRR7211988", "depressus", "droebachiensis", "fragilis", "franciscanus", "intermedius", "nudus", "pallidus", "pulcherrimus_DRR107784", "pulcherrimus_SRR5767283", "purpuratus_SRR6281818", "purpuratus_SRR7211988"]
	writer.writerow(header)

	for key,value in passed_genes_dict.items():
		row = [key] + value[0] + value[1]
		writer.writerow(row)

	csv_file.close()

def main():
	bed_file_list = get_zipped_bed_file_list()
	
	initialize_gene_dict()

	for regions_file, thresholds_file in bed_file_list:
		fill_gene_dict(regions_file, thresholds_file)

	filter_gene_dict()

	write_genes_passed_filter_bed()

	write_all_gene_dict_csv()

	write_passed_genes_dict_csv()

if __name__ == "__main__":
	main()
