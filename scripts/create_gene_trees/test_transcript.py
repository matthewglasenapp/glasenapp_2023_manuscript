import os
import gzip
from itertools import islice
import csv
from joblib import Parallel, delayed
import multiprocessing
import subprocess

num_cores = multiprocessing.cpu_count()

# S. purpuratus gff3 file
gff_file = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Bed file containing a record for each protein coding gene in the S. purpuratus assembly. See the ncbi/ directory for scripts to generate this file
protein_coding_genes_bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_genes/protein_coding_genes.bed"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
# Strongylocentrotidae Subset
subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Average coverage of S. purpuratus genes for each sample
mean_coverage_spur5_genes = {
"depressus_SRR5767284": 22,
"droebachiensis_SRR5767286" : 26.46,
"fragilis_SRR5767279": 33.73,
"franciscanus_SRR5767282": 22.47,
"intermedius_SRR5767280": 30.95,
"lividus_ERS2351987": 6.22,
"nudus_SRR5767281": 24.1,
"pallidus_SRR5767285": 12.01,
"pulcherrimus_DRR107784": 77.65,
"pulcherrimus_SRR5767283": 28.15,
"purpuratus_SRR6281818": 51.56,
"purpuratus_SRR7211988": 93.9,
"variegatus_SRR7207203": 8.4
}

# Average coverage of S. purpuratus exons for each sample
mean_coverage_spur5_exons = {
"depressus_SRR5767284": 47.5,
"droebachiensis_SRR5767286" : 41.5,
"fragilis_SRR5767279": 46.8,
"franciscanus_SRR5767282": 33.8,
"intermedius_SRR5767280": 44.2,
"lividus_ERS2351987": 12,
"nudus_SRR5767281": 40.5,
"pallidus_SRR5767285": 15,
"pulcherrimus_DRR107784": 108.47,
"pulcherrimus_SRR5767283": 44.3,
"purpuratus_SRR6281818": 55.81,
"purpuratus_SRR7211988": 100.3,
"variegatus_SRR7207203": 8.4
}

#bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/unique_exons.bed"
bed_file = "/hb/scratch/mglasena/mrna_cov_2/pallidus_SRR5767285.regions.bed.gz"

#bed_file_dir = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_genes/"
#bed_file_dir = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/"
bed_file_dir = "/hb/scratch/mglasena/mrna_cov_2/"

min_cov_threshold = 10

# Consider not filtering by prop_1x_threshold, because this value is heavily determined by what proportion of the gene is actually CDS. 
# Enter as a proportion
prop_1x_threshold = 0.9

prop_10x_threshold = 0.0

subset_mean_coverage_spur5_exons = dict()

# {"rna": ("parent_gene", [average coverage depth for each species], [prop_1x], [prop_10x], [prop_20x])}
rna_dict = dict()

passed_genes_dict = dict()

# For filtering by gap part of code
#required_gap = 100000
required_gap = 50000
# Scaffold dictionay in format of "Scaffold":[[scaffold,gene,start,stop]]
scaffold_dict = {}
# Initialize list for genes passing gap filter
passed_genes = []

# Dictionary linking mRNA to parent genes
mrna_gene_dict = dict()

filtered_mrna_gene_dict = dict()

# For vcf2fasta part of the code 
vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"
vcf_file = "/hb/scratch/mglasena/data/genotypes/strongylocentrotidae/3bp_filtered_genotype_calls_pf.g.vcf.gz"
feature = "gene"

sample_names = {
#'(4': "(lividus",
'QB3KMK013': 'fragilis',
'QB3KMK011': 'nudus',
'QB3KMK010': 'franciscanus',
'QB3KMK015': 'depressus',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
#'QB3KMK016': 'pulcherrimus_SRR5767283',
'QB3KMK016': 'pulcherrimus',
'QB3KMK012': 'intermedius',
#'SPUR.00': 'purpuratus_SRR7211988',
'SPUR.00': 'purpuratus',
#'SAMD00098133': 'pulcherrimus_DRR107784',
#'S.purpuratus_1': 'purpuratus_SRR6281818',
#'LVAR.00': 'variegatus'
}

# nw_utils directory
nw_utils = "/hb/groups/pogson_group/dissertation/software/newick_utils/src/"

def subset_coverage_dict():
	try:
		for sample in subset_sample_list:
			subset_mean_coverage_spur5_exons[sample] = mean_coverage_spur5_exons[sample]
	
	except:
		for key in mean_coverage_spur5_exons:
			subset_mean_coverage_spur5_exons[key] = mean_coverage_spur5_exons[key]

def get_zipped_bed_file_list():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	return list(file_list)

# Create dictionary in the format of 
def initialize_gene_dict():
	with gzip.open(bed_file,"rt") as f:
		for line in f:
			gene = line.split("\t")[0]
			rna_dict[gene] = [],[],[],[]

def fill_gene_dict(regions_file, thresholds_file):
	with gzip.open(regions_file, "rt") as file1, gzip.open(thresholds_file, "rt") as file2:
		for line_file1, line_file2 in zip(file1, islice(file2, 0, None)):
			gene_id = line_file1.split("\t")[-2]
			mean_depth = float(line_file1.split("\t")[-1].strip())

			total_base_count = int(line_file2.split("\t")[-4])
			
			one_x_count = int(line_file2.split("\t")[-3])
			ten_x_count = int(line_file2.split("\t")[-2])
			twenty_x_count = int(line_file2.split("\t")[-1])
			
			try:
				prop_1x = float(one_x_count / total_base_count)
				prop_10x = float(ten_x_count / total_base_count)
				prop_20x = float(twenty_x_count / total_base_count)
			except ZeroDivisionError:
				prop_1x = 0.0
				prop_10x = 0.0
				prop_20x = 0.0

			rna_dict[gene_id][0].append(mean_depth)
			rna_dict[gene_id][1].append(prop_1x)
			rna_dict[gene_id][2].append(prop_10x)
			rna_dict[gene_id][3].append(prop_20x)

def write_all_rna_dict_csv():
	csv_file = open("all_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	try:
		for i in range(4):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(4):
			for key in mean_coverage_spur5_genes.keys():
				header.append(key)
	writer.writerow(header)

	for key,value in rna_dict.items():
		row = [key] + value[0] + value[1] + value[2] + value[3]
		writer.writerow(row)

	csv_file.close()

def filter_gene_dict():
	for key, value in rna_dict.items():
		mean_depth_lst = [item for item in value[0]]
		one_x_lst = [item for item in value[1]]
		ten_x_lst = [item for item in value[2]]
		twenty_x_lst = [item for item in value[3]]
		
		if min(mean_depth_lst) >= min_cov_threshold and min(one_x_lst) >= prop_1x_threshold and min(ten_x_lst) >= prop_10x_threshold:
			test_var = True
			for counter, sample in enumerate(subset_mean_coverage_spur5_exons):
				if mean_depth_lst[counter] >= (subset_mean_coverage_spur5_exons[sample] * 2):
					test_var = False
				
			if test_var == True:
				passed_genes_dict[key] = value

def write_passed_genes():
	get_mrna_from_gff = '''awk '$3 == "mRNA"' {} > mrna_records.txt'''.format(gff_file)
	os.system(get_mrna_from_gff)
	
	# Creates mapping of rna: parent_gene for mRNAs in passed_genes_dict 
	with open("mrna_records.txt", "r") as f:
		for line in f:
			scaffold = line.split("\t")[0]
			start = line.split("\t")[3]
			stop = line.split("\t")[4]
			rna = line.split("\t")[8].split(";")[0].split("rna-")[1]
			parent_gene = line.split("\t")[8].split(";")[1].split("gene-")[1]
			
			if rna in passed_genes_dict.keys() and scaffold != "NC_001453.1":
				mrna_gene_dict[rna] = parent_gene

				if scaffold_dict.get(scaffold):
					scaffold_dict[scaffold].append([scaffold,rna,start,stop])
				else:
					scaffold_dict[scaffold] = [[scaffold,rna,start,stop]]

def check_proximity():
	gene_counter = 0 
	for gene_list in scaffold_dict.values():
		for gene in gene_list:
			try:
				current_stop = gene_list[gene_counter][3]

			except IndexError:
				break
				
			try:
				next_start = gene_list[gene_counter + 1][2]
			
			except IndexError:
				gene_counter = 0
				break
			
			while (int(next_start) - int(current_stop)) < required_gap:
				current_gene = gene_list[gene_counter]
				failed_gene = gene_list.pop(gene_counter+1)

				try:
					next_start = gene_list[gene_counter+1][2]
				except IndexError:
					break
			
			gene_counter += 1

def get_passed_genes_list():
	for gene_list in scaffold_dict.values():
		for gene in gene_list:
			passed_genes.append(gene[1])

	for key,value in mrna_gene_dict.items():
		if key in passed_genes:
			filtered_mrna_gene_dict[key] = value

def write_new_bed_file():
	with open(protein_coding_genes_bed_file,"r") as f:
		with open("unlinked_loci.bed","a") as f2:
			for line in f:
				if line.split("\t")[3].split("gene-")[1] in filtered_mrna_gene_dict.values():
					f2.write(line)

def write_passed_gene_dict_csv():

	csv_file = open("passed_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	try:
		for i in range(4):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(4):
			for key in mean_coverage_spur5_genes.keys():
				header.append(key)
	writer.writerow(header)

	for key,value in passed_genes_dict.items():
		if key in filtered_mrna_gene_dict.keys():
			row = [key] + value[0] + value[1] + value[2] + value[3]
			writer.writerow(row)

	csv_file.close()

def get_gene_ids():
	split_columns = "awk '{ print $10 }' unlinked_loci.bed > gene_list"
	os.system(split_columns)

	with open("gene_list","r") as f:
		with open("gene_ids","a") as f2:
			gene_list = f.read().splitlines()
			for gene in gene_list:
				identifier = gene.split(";")[1]
				# For exons
				#identifier = gene.split("\t")[3]
				f2.write(identifier + "\n")

	gene_ids = open("gene_ids", "r").read().splitlines()
	return gene_ids

def make_sco_gff(gene):
	command = "grep {} {} > {}.record".format(gene, gff_file, gene)
	os.system(command)

def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff sco_gff.gff --feat {}".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)
	
def run_iqtree():
	#run_iqtree = "iqtree -S vcf2fasta_gene/ -m MFP --prefix loci -T AUTO"
	run_iqtree = "iqtree -S vcf2fasta_gene/ -m GTR -o QB3KMK016 --prefix loci -T AUTO -B 1000 --boot-trees"
	os.system(run_iqtree)

def subset_boot_file():
	line_counter = 0 

	with open("loci.ufboot_subset","a") as f:
		for line in open("loci.ufboot","r"):
			if ((line_counter + 1000) % 1000) == 0:
				count = 0
			if count <= 99:
				f.write(line)
				count += 1 
			line_counter += 1

def edit_tree_files(input_file, output_file):
	with open(input_file, "r") as f:
		tree_list = f.read().splitlines()
	
	with open(output_file,"a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

def clean_gene_trees(input_file, output_file):
	clean = "{}nw_topology -I {} | {}nw_order -c d - > {}".format(nw_utils, input_file, nw_utils, output_file)
	os.system(clean)

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


	write_all_rna_dict_csv()

	filter_gene_dict()

	write_passed_genes()
	
	check_proximity()
	get_passed_genes_list()
	
	write_new_bed_file()
	write_passed_gene_dict_csv()

	#gene_ids = get_gene_ids()

	#Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	#os.system("cat *.record > sco_gff.gff")
	#os.system("rm *.record")

	#run_vcf2fasta()
	#replace_missing_genotype_char()
	#run_iqtree()
	#subset_boot_file()
	#edit_tree_files("loci.treefile","single_locus_trees.nwk")
	#edit_tree_files("loci.ufboot_subset", "single_locus_trees_boot_subset.nwk")
	#clean_gene_trees("single_locus_trees.nwk", "clean_trees.nwk")

if __name__ == "__main__":
	main()