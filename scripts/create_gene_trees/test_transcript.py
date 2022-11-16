import os
import gzip
from itertools import islice
import csv
from joblib import Parallel, delayed
import multiprocessing
import subprocess

num_cores = multiprocessing.cpu_count()

# nw_utils directory
nw_utils = "/hb/groups/pogson_group/dissertation/software/newick_utils/src/"

# Path to vcf2fasta.py
vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"

# Feature of gff file that vcf2fasta.py will build alignments for
feature = "CDS"

# Path to S. purpuratus reference genome
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# S. purpuratus gff3 file
gff_file = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Path to filtered multisample vcf file
vcf_file = "/hb/scratch/mglasena/data/genotypes/strongylocentrotidae/3bp_filtered_genotype_calls_pf.g.vcf.gz"

# Bed file containing a record for each protein coding gene in the S. purpuratus assembly. See the ncbi/ directory for scripts to generate this file
protein_coding_genes_bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_genes/protein_coding_genes.bed"

# Directory containing the output files from mrna.py
bed_file_dir = "/hb/scratch/mglasena/mrna_cov_2/"

# Sample bed file from the output of mrna.py
bed_file = "/hb/scratch/mglasena/mrna_cov_2/pallidus_SRR5767285.regions.bed.gz"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
# Strongylocentrotidae Subset
subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Specify thresholds for filtering. 
min_cov_threshold = 10

prop_1x_threshold = 0.9

prop_10x_threshold = 0.0

# Required gap between two adjacent loci in base pairs 
required_gap = 50000

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

# Mapping of DNA sample names to species names 
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

# Initialize new dictionary containing the average coverage of S. purpuratus exons for each sample in subset_sample_list variable
subset_mean_coverage_spur5_exons = dict()

# Initialize dictionary in format of {"rna": [average coverage depth for each species], [prop_1x], [prop_10x], [prop_20x])}
rna_dict = dict()

# Initialize dictionary for mRNAs passing initial filter by coverage depth 
passed_rna_dict = dict()

# Initialize dictionary mapping mRNAs passing intial filter by coverage depth to the name of their parent gene
mrna_gene_dict = dict()

# Initialize dictionay in format of "Scaffold":[[scaffold,rna,start,stop]] for mRNAs passing initial filter. 
# This dictionary will be used to filter out mRNAs that are too close to one another
scaffold_dict = {}

# Initialize final list for genes passing both the coverage depth and gap filters
passed_rnas = []

# Initialize dictionary mapping mRNA names to parent gene name for mRNAs that passed all filters
filtered_mrna_gene_dict = dict()

# Populate subset_mean_coverage_spur5_exons variable with the samples included in the subset_sample_list variable
def subset_coverage_dict():
	try:
		for sample in subset_sample_list:
			subset_mean_coverage_spur5_exons[sample] = mean_coverage_spur5_exons[sample]
	
	except:
		for key in mean_coverage_spur5_exons:
			subset_mean_coverage_spur5_exons[key] = mean_coverage_spur5_exons[key]

# Get zipped list of regions and thresholds files for each species 
def get_zipped_bed_file_list():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	os.system("rm regions_files")
	os.system("rm thresholds_files")

	return list(file_list)

# Populate keys and create structure for rna_dict in format of {"rna": [average coverage depth for each species], [prop_1x], [prop_10x], [prop_20x])}
def initialize_rna_dict():
	with gzip.open(bed_file,"rt") as f:
		for line in f:
			gene = line.split("\t")[0]
			rna_dict[gene] = [],[],[],[]

# Populate rna dict with coverage metrics for each species. No filtering involved in this step 
def fill_rna_dict(regions_file, thresholds_file):
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

# Write rna_dict to csv file 
def write_all_rna_dict_csv():
	csv_file = open("all_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	
	# Create header row for csv file
	try:
		for i in range(4):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(4):
			for key in mean_coverage_spur5_genes.keys():
				header.append(key)
	
	writer.writerow(header)

	# Write row for each mRNA containing mean depth, prop1x, prop10x, and prop20x for each species/sample
	counter = 0
	for key,value in rna_dict.items():
		row = [key] + value[0] + value[1] + value[2] + value[3]
		counter += 1
		writer.writerow(row)

	print("{} mRNA records written to all_rna.csv".format(counter))
	csv_file.close()

# Filter rna_dict by coverage metrics. Add mRNA's passing filter to passed_rna_dict
def filter_rna_dict():
	for key, value in rna_dict.items():
		mean_depth_lst = [item for item in value[0]]
		one_x_lst = [item for item in value[1]]
		ten_x_lst = [item for item in value[2]]
		twenty_x_lst = [item for item in value[3]]
		
		if min(mean_depth_lst) >= min_cov_threshold and min(one_x_lst) >= prop_1x_threshold and min(ten_x_lst) >= prop_10x_threshold:
			test_var = True
			
			# Filter any species who's coverage for that mRNA is double or greater the average coverage for exons
			for counter, sample in enumerate(subset_mean_coverage_spur5_exons):
				if mean_depth_lst[counter] >= (subset_mean_coverage_spur5_exons[sample] * 2):
					test_var = False
				
			if test_var == True:
				passed_rna_dict[key] = value

	print("{} mRNA records passed intial coverage depth filters".format(len(passed_rna_dict)))

# Get file of mRNA records from gff3 file
def get_mrna_gff():
	get_mrna_from_gff = '''awk '$3 == "mRNA"' {} > mrna_records.txt'''.format(gff_file)
	os.system(get_mrna_from_gff)
	
# Create dictionary for mRNAs passing initial filter in the format of "Scaffold":[[scaffold,rna,start,stop]]. This dictionary will be used for filtering by distance on chromosome
def create_scaffold_dict():
	
	# Creates mapping of rna: parent_gene for mRNAs in passed_rna_dict 
	with open("mrna_records.txt", "r") as f:
		for line in f:
			scaffold = line.split("\t")[0]
			start = line.split("\t")[3]
			stop = line.split("\t")[4]
			rna = line.split("\t")[8].split(";")[0].split("rna-")[1]
			parent_gene = line.split("\t")[8].split(";")[1].split("gene-")[1]
			
			mt_rna_counter = 0 
			if rna in passed_rna_dict.keys():
				
				if scaffold != "NC_001453.1":
					mrna_gene_dict[rna] = parent_gene
				
				else:
					mt_rna_counter += 1

				if scaffold_dict.get(scaffold):
					scaffold_dict[scaffold].append([scaffold,rna,start,stop])
				else:
					scaffold_dict[scaffold] = [[scaffold,rna,start,stop]]

	record_counter = 0
	for value in scaffold_dict.values():
		for item in value:
			record_counter += 1

	print(scaffold_dict)
	os.system("rm mrna_records.txt")
	print("{} mitochondrial mRNAs were removed".format(mt_rna_counter))
	print("{} mRNAs added to scaffold_dict".format(record_counter))

# Filter mRNA records in scaffold_dict so that all remaining mRNAs are >= required_gap apart from each other
def check_proximity():
	rna_counter = 0
	filter_counter = 0
	
	for rna_list in scaffold_dict.values():
		for rna in rna_list:
			try:
				current_stop = rna_list[rna_counter][3]

			except IndexError:
				break
				
			try:
				next_start = rna_list[rna_counter + 1][2]
			
			except IndexError:
				rna_counter = 0
				break
			
			while (int(next_start) - int(current_stop)) < required_gap:
				current_rna = rna_list[rna_counter]
				rna_list.pop(rna_counter + 1)
				filter_counter += 1

				try:
					next_start = rna_list[rna_counter + 1][2]
				
				except IndexError:
					break
			
			rna_counter += 1

	record_counter = 0
	for value in scaffold_dict.values():
		for item in value:
			record_counter += 1

	print("{} mRNAs removed to satisfy gap filter".format(filter_counter))
	print("{} mRNAs passed both coverage depth and gap filters".format(record_counter))

# Get list of passed mRNAs and create new filtered_mrna_gene_dict that only includes mRNAs that passed the gap filter. 
def get_passed_rnas():
	for rna_list in scaffold_dict.values():
		for rna in rna_list:
			passed_rnas.append(rna[1])

	for key,value in mrna_gene_dict.items():
		if key in passed_rnas:
			filtered_mrna_gene_dict[key] = value

	print("{} mRNAs in passed_rnas list".format(len(passed_rnas)))
	print("{} mRNAs in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))

	with open("test_list","a") as f:
		for record in passed_rnas:
			f.write(record + "\n")


# Write new bed file ("unlinked_loci.bed") of parent genes with an mRNA transcript that passed all filters
def write_new_bed_file():
	records_written = 0
	records_written_lst = []
	
	with open(protein_coding_genes_bed_file,"r") as f:
		with open("unlinked_loci.bed","a") as f2:
			for line in f:
				if line.split("\t")[3].split("gene-")[1] in filtered_mrna_gene_dict.values():
					records_written_lst.append(line.split("\t")[3].split("gene-")[1])
					records_written += 1
					f2.write(line)

	print("{} records written to unlinked_loci.bed".format(records_written))

# Write csv file of RNAs passing all filters and their coverage metrics
def write_passed_rna_dict_csv():
	csv_file = open("passed_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	records_written = 0 
	
	try:
		for i in range(4):
			for sample in subset_sample_list:
				header.append(sample)
	except NameError:
		for i in range(4):
			for key in mean_coverage_spur5_genes.keys():
				header.append(key)
	
	writer.writerow(header)

	for key,value in passed_rna_dict.items():
		if key in filtered_mrna_gene_dict.keys():
			row = [key] + value[0] + value[1] + value[2] + value[3]
			records_written += 1
			writer.writerow(row)

	csv_file.close()

	print("{} records written to passed_rna.csv".format(records_written))

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
	
	initialize_rna_dict()

	for regions_file, thresholds_file in bed_file_list:
		try:
			for sample in subset_sample_list:
				if sample in regions_file and sample in thresholds_file:
					fill_rna_dict(regions_file, thresholds_file)
		
		except NameError:
			fill_rna_dict(regions_file, thresholds_file)

	print("There were {} mRNA records pre-filter".format(len(rna_dict)))

	write_all_rna_dict_csv()
	filter_rna_dict()
	get_mrna_gff()
	create_scaffold_dict()
	check_proximity()
	get_passed_rnas()
	write_new_bed_file()
	write_passed_rna_dict_csv()

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