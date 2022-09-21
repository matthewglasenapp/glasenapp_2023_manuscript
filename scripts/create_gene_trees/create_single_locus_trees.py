# Two gff_file objects!

import os
import gzip
from itertools import islice
import csv
from joblib import Parallel, delayed
import multiprocessing
import subprocess

num_cores = multiprocessing.cpu_count()

gff_file = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
#subset_sample_list = ["droebachiensis_SRR5767286", "fragilis_SRR5767279", "intermedius_SRR5767280", "pallidus_SRR5767285"]
#subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_DRR107784', 'pulcherrimus_SRR5767283', 'purpuratus_SRR6281818', 'purpuratus_SRR7211988']

mean_coverage_spur5_exons = {
"depressus_SRR5767284": 47.45,
"droebachiensis_SRR5767286" : 41.46,
"fragilis_SRR5767279": 46.75,
"franciscanus_SRR5767282": 33.77,
"intermedius_SRR5767280": 44.21,
"lividus_ERS2351987": 11.99,
"nudus_SRR5767281": 40.53,
"pallidus_SRR5767285": 15,
"pulcherrimus_DRR107784": 108.47,
"pulcherrimus_SRR5767283": 44.34,
"purpuratus_SRR6281818": 55.81,
"purpuratus_SRR7211988": 100.31
}

bed_file = "/hb/scratch/mglasena/run_mosdepth/protein_coding_genes.bed"

bed_file_dir = "/hb/scratch/mglasena/run_mosdepth/"

min_cov_threshold = 5

# Consdier not filtering by prop_1x_threshold, because this value is heavily determined by what proportion of the gene is actually CDS. 
prop_1x_threshold = 0.0

subset_mean_coverage_spur5_exons = dict()

gene_dict = dict()

passed_genes_dict = dict()

# For vcf2fasta part of the code 
vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"
vcf_file = "/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls.g.vcf.gz"
gff_file = "/hb/scratch/mglasena/test/sco_gff.gff"
feature = "gene"

sample_names = {
'QB3KMK013': 'fragilis',
#'QB3KMK011': 'nudus',
'QB3KMK010': 'franciscanus',
#'QB3KMK015': 'depressus',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
#'S.purpuratus#1': 'purpuratus_SRR6281818',
'QB3KMK016': 'pulcherrimus_SRR5767283',
'QB3KMK012': 'intermedius',
'SPUR.00': 'purpuratus_SRR7211988',
#'SAMD00098133': 'pulcherrimus_DRR107784'
}

def subset_coverage_dict():
	try:
		for sample in subset_sample_list:
			subset_mean_coverage_spur5_exons[sample] = mean_coverage_spur5_exons[sample]
	
	except:
		for key in mean_coverage_spur5_exons:
			subset_mean_coverage_spur5_exons[key] = mean_coverage_spur5_exons[key]

def get_zipped_bed_file_list():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* | grep -v 'csi' > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	return list(file_list)

# Create dictionary in the format of {"gene_id": [average coverage depth for each species], [proportion of bases covered for each species]}
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
		mean_depth_lst = [item for item in value[0]]
		one_x_lst = [item for item in value[1]]
		
		if min(mean_depth_lst) >= min_cov_threshold and min(one_x_lst) >= prop_1x_threshold:
			test_var = True
			for counter, sample in enumerate(subset_mean_coverage_spur5_exons):
				if mean_depth_lst[counter] >= (subset_mean_coverage_spur5_exons[sample] * 2):
					test_var = False
				
			if test_var == True:
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
	header = ["gene"]
	try:
		for i in range(2):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(2):
			for key in mean_coverage_spur5_exons.keys():
				header.append(key)
	writer.writerow(header)

	for key,value in gene_dict.items():
		row = [key] + value[0] + value[1]
		writer.writerow(row)

	csv_file.close()

def write_passed_genes_dict_csv():
	csv_file = open("passed_genes.csv","w")
	writer = csv.writer(csv_file)	
	header = ["gene"]
	try:
		for i in range(2):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(2):
			for key in mean_coverage_spur5_exons.keys():
				header.append(key)
	writer.writerow(header)

	for key,value in passed_genes_dict.items():
		row = [key] + value[0] + value[1]
		writer.writerow(row)

	csv_file.close()

def get_gene_ids():
	split_columns = "awk '{ print $10 }' genes_pass_filter.bed > gene_list"
	os.system(split_columns)

	with open("gene_list","r") as f:
		with open("gene_ids","a") as f2:
			gene_list = f.read().splitlines()
			for gene in gene_list:
				identifier = gene.split(";")[1]
				f2.write(identifier + "\n")

	gene_ids = open("gene_ids", "r").read().splitlines()
	return gene_ids

def make_sco_gff(gene):
	command = "grep {} {} > {}.txt".format(gene, gff_file, gene)
	os.system(command)

def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff sco_gff --feat {}".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)
	
def run_iqtree():
	run_iqtree = "iqtree -S vcf2fasta_gene/ -m MFP --prefix loci -T AUTO"
	os.system(run_iqtree)

def edit_tree_files():
	with open("loci.treefile", "r") as f:
		tree_list = f.read().splitlines()
	
	with open("single_locus_trees.nwk","a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

def main():
	subset_coverage_dict()

	bed_file_list = get_zipped_bed_file_list()
	
	initialize_gene_dict()

	for regions_file, thresholds_file in bed_file_list:
		try:
			for sample in subset_sample_list:
				if sample in regions_file and sample in thresholds_file:
					fill_gene_dict(regions_file, thresholds_file)
		
		except NameError:
			fill_gene_dict(regions_file, thresholds_file)

	filter_gene_dict()

	write_genes_passed_filter_bed()

	write_all_gene_dict_csv()

	write_passed_genes_dict_csv()

	gene_ids = get_gene_ids()

	Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	os.system("cat *.txt > sco_gff")
	os.system("rm *.txt")

	run_vcf2fasta()
	replace_missing_genotype_char()
	run_iqtree()
	edit_tree_files()

if __name__ == "__main__":
	main()