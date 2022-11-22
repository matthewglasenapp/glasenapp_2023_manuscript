import os
import gzip
from itertools import islice
import csv
from joblib import Parallel, delayed
import multiprocessing
import subprocess
import statistics
import multiprocessing

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
vcf_file = "/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls_pf.g.vcf.gz"

# Bed file containing a record for each protein coding gene in the S. purpuratus assembly. See the ncbi/ directory for scripts to generate this file
protein_coding_genes_bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_genes/protein_coding_genes.bed"

# Directory containing the output files from mrna.py
bed_file_dir = "/hb/scratch/mglasena/mrna_cov/"

# Sample bed file from the output of mrna.py
bed_file = "/hb/scratch/mglasena/mrna_cov/pallidus_SRR5767285.regions.bed.gz"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
# Strongylocentrotidae Subset
#subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Franciscanus outgroup
subset_sample_list = ['droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Specify thresholds for filtering. 
min_cov_threshold = 15

prop_1x_threshold = 0.9

prop_10x_threshold = 0.9

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

# Initialize dictionary mapping CDS names to parent mRNA names for CDS of genes that had an mRNA pass previous filters. 
cds_parent_rna_dict = dict()

# List of redundant isoforms to delete following vcf2fasta
records_to_delete = []

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

	os.system("rm mrna_records.txt")
	print("{} mitochondrial mRNAs were removed".format(mt_rna_counter))
	print("{} mRNAs added to scaffold_dict".format(record_counter))

# Filter mRNA records in scaffold_dict so that all remaining mRNAs are >= required_gap apart from each other
def check_proximity():
	filter_counter = 0
	failed_rna_lst = []
	
	for rna_list in scaffold_dict.values():
		first = True
		last_stop = 0
		rna_counter = 0

		while rna_counter < len(rna_list):
			(chromosome, gene, start, stop) = rna_list[rna_counter]

			if first:
				last_stop = int(stop)
				first = False
				rna_counter += 1
				continue
			
			elif int(start) < (last_stop + required_gap):

				failed_rna = rna_list.pop(rna_counter)
				#failed_rna_lst.append(failed_rna)
				filter_counter +=1 
			
			else:
				last_stop = int(stop)
				rna_counter +=1

	#print(failed_rna_lst)
	
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

	with open("passed_rnas.txt","a") as f:
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



# Get list of parent gene identifiers for those genes that passed all filters. Example: Dbxref=GeneID:582406
def get_gene_ids():
	get_info_column = "awk '{ print $10 }' unlinked_loci.bed > gene_list"
	os.system(get_info_column)

	with open("gene_list","r") as f, open("gene_ids","a") as f2:
		gene_list = f.read().splitlines()
		for gene in gene_list:
			identifier = gene.split(";")[1]
			f2.write(identifier + "\n")

	os.system("rm gene_list")
	gene_ids = open("gene_ids", "r").read().splitlines()
	os.system("rm gene_ids")
	return gene_ids

# Using S. purpuratus gff3 file, make gff file for each gene that passed previous filters
def make_sco_gff(gene):
	command = "grep {} {} > single_gene_gff_records/{}.record".format(gene, gff_file, gene)
	os.system(command)

# Use vcf2fasta to create fasta alignments for all mRNAs passing filter
def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff sco_gff.gff --feat {} --blend".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

# Vcf2fasta makes alignments for all isoforms of each gene. This function identifies the unwanted isoform alignment files and adds them to the records_to_delete list.
def identify_redundant_isoforms():
	get_cds_gff = '''awk '$3 == "CDS"' sco_gff.gff > cds.gff'''
	os.system(get_cds_gff)

	with open("cds.gff", "r") as f:
		records = f.read().splitlines()

		for record in records:
			cds_name = record.split("\t")[8].split(";")[0].split("cds-")[1]
			parent_rna_name = record.split("\t")[8].split(";")[1].split("rna-")[1]

			if cds_name in cds_parent_rna_dict:
				if cds_parent_rna_dict[cds_name] == parent_rna_name:
					continue
			
				else:
					print(cds_name)
					print(parent_rna_name)
					break

			else:
				cds_parent_rna_dict[cds_name] = parent_rna_name

	with open("passed_rnas.txt") as f2:
		passed_rnas_lst = f2.read().splitlines()

		for key,value in cds_parent_rna_dict.items():
			if value in passed_rnas_lst:
				continue

			else:
				records_to_delete.append(key)

	os.system("rm sco_gff.gff")
	os.system("rm cds.gff")

# Delete redundant isoform alignment files 
def delete_redundant_isoforms():
	for record in records_to_delete:
		delete = "rm vcf2fasta_CDS/{}.fas".format(record)
		os.system(delete)

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_CDS/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_CDS/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

# Run iqturee on the fasta files of mRNAs that passed all filters
def identify_no_variant_no_parsimony():
	#run_iqtree = "iqtree2 -S vcf2fasta_CDS/ -m MFP --prefix loci -T AUTO -B 1000 --boot-trees"
	#run_iqtree = "iqtree2 -S vcf2fasta_CDS/ -m MFP --prefix loci -T 8"
	run_iqtree = "iqtree2 -S vcf2fasta_CDS/ -m TESTONLY --prefix loci -T 8"
	os.system(run_iqtree)

def remove_no_variant_no_parsimony():
	no_parsimony_lst = []
	no_variant_lst = []
	
	get_no_parsimony = '''cat loci.log | grep "No parsimony" > no_parsimony.txt'''
	get_no_variant = '''cat loci.log | grep "No variant" > no_variant.txt'''

	os.system(get_no_parsimony)
	os.system(get_no_variant)

	with open("no_parsimony.txt", "r") as f, open("no_variant.txt", "r") as f2:
		for line in f:
			no_parsimony_lst.append(line.split(" ")[6].strip())

		for line in f2:
			no_variant_lst.append(line.split(" ")[6].strip())

	os.mkdir('no_variant')
	os.mkdir('no_parsimony')

	for file in no_parsimony_lst:
		move = "mv vcf2fasta_CDS/{} no_parsimony/".format(file)
		os.system(move)

	for file in no_variant_lst:
		move = "mv vcf2fasta_CDS/{} no_variant/".format(file)
		os.system(move)

	os.system("rm no_parsimony.txt")
	os.system("rm no_variant.txt")
	os.system("rm *loci*")

def run_iqtree(fasta_file):
	run_iqtree = "iqtree2 -s vcf2fasta_CDS/{} -m MFP -b 100 --boot-trees -T 2".format(fasta_file)
	os.system(run_iqtree)

def clean_up_iqtree_files():
	delete_gz = 'find ./vcf2fasta_CDS/ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./vcf2fasta_CDS/ -type f -name "*.bionj" -delete'
	delete_contree = 'find ./vcf2fasta_CDS/ -type f -name "*.contree" -delete'
	delete_iqtree = 'find ./vcf2fasta_CDS/ -type f -name "*.iqtree" -delete'
	delete_log = 'find ./vcf2fasta_CDS/ -type f -name "*.log" -delete'
	delete_mldist = 'find ./vcf2fasta_CDS/ -type f -name "*.mldist" -delete'
	os.system(delete_gz)
	os.system(delete_bionj)
	os.system(delete_contree)
	os.system(delete_iqtree)
	os.system(delete_log)
	os.system(delete_mldist)

	cat_treefiles = 'find ./vcf2fasta_CDS/ -type f -name "*.treefile" -exec cat {} \\; > loci.treefile'
	cat_boottrees = 'find ./vcf2fasta_CDS/ -type f -name "*.boottrees" -exec cat {} \\; > loci.boottrees'
	
	os.system(cat_treefiles)
	os.system(cat_boottrees)

	delete_treefile = 'find ./vcf2fasta_CDS/ -type f -name "*.treefile" -delete'
	delete_boottrees = 'find ./vcf2fasta_CDS/ -type f -name "*.boottrees" -delete'
	
	os.system(delete_treefile)
	os.system(delete_boottrees)

# Subset the --boot-trees file produced by iqtree. Currently, -b does not work with -S. -B requires >= 1000. I only want 100 boot trees per locus though. 
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

# Edit tree files to show species names instead of sample names. 
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

# Clean trees to remove branch lengths and bootstrap support values. Arrange identical topolgoies to have the same textual representation
def clean_gene_trees(input_file, output_file):
	clean = "{}nw_topology -I {} | {}nw_order -c d - > {}".format(nw_utils, input_file, nw_utils, output_file)
	os.system(clean)

def main():
	#subset_coverage_dict()

	#bed_file_list = get_zipped_bed_file_list()
	
	#initialize_rna_dict()

	#for regions_file, thresholds_file in bed_file_list:
		#try:
			#for sample in subset_sample_list:
				#if sample in regions_file and sample in thresholds_file:
					#fill_rna_dict(regions_file, thresholds_file)
		
		#except NameError:
			#fill_rna_dict(regions_file, thresholds_file)

	#print("There were {} mRNA records pre-filter".format(len(rna_dict)))

	#write_all_rna_dict_csv()
	#filter_rna_dict()
	#get_mrna_gff()
	#create_scaffold_dict()
	#check_proximity()
	#get_passed_rnas()
	#write_new_bed_file()
	#write_passed_rna_dict_csv()

	#gene_ids = get_gene_ids()

	#os.system("mkdir single_gene_gff_records/")
	#Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	
	# Concatenate all single gene gff records into "sco_gff.gff" file
	#os.system('find ./single_gene_gff_records/ -type f -name "*.record" -exec cat {} \\; > sco_gff.gff')
	
	# Delete the single gene records
	#os.system('find ./single_gene_gff_records/ -type f -name "*.record" -delete')
	#os.system('rmdir single_gene_gff_records/')

	#run_vcf2fasta()

	#identify_redundant_isoforms()
	#delete_redundant_isoforms()

	#replace_missing_genotype_char()
	#identify_no_variant_no_parsimony()
	#remove_no_variant_no_parsimony()
	
	#fasta_file_list = os.listdir("vcf2fasta_CDS")[0:872]
	#fasta_file_list = os.listdir("vcf2fasta_CDS")[872:1744]
	#fasta_file_list = os.listdir("vcf2fasta_CDS")[1744:2616]

	fasta_files = ['XP_030833507.1.fas', 'XP_030853234.1.fas', 'XP_030835274.1.fas', 'XP_030840120.1.fas', 'XP_030838814.1.fas', 'XP_793152.2.fas', 'XP_030835872.1.fas', 'XP_030832482.1.fas', 'XP_030845105.1.fas', 'XP_030853413.1.fas', 'XP_783931.2.fas', 'XP_030837783.1.fas', 'XP_030846089.1.fas', 'XP_030836903.1.fas', 'XP_786209.3.fas', 'XP_030830175.1.fas', 'XP_003728452.2.fas', 'XP_030832258.1.fas', 'XP_011668113.1.fas', 'XP_030837711.1.fas', 'XP_030840768.1.fas', 'XP_011681581.1.fas', 'XP_030839903.1.fas', 'XP_011683286.2.fas', 'XP_030832113.1.fas', 'XP_011676429.1.fas', 'XP_011677085.2.fas', 'XP_030844235.1.fas', 'XP_011683752.1.fas', 'XP_030854389.1.fas', 'XP_011674607.1.fas', 'XP_003727052.2.fas', 'XP_784087.1.fas', 'XP_011675281.1.fas', 'XP_030828856.1.fas', 'XP_030856312.1.fas', 'XP_030846855.1.fas', 'XP_030839892.1.fas', 'XP_030845006.1.fas', 'XP_030833404.1.fas', 'XP_030852849.1.fas', 'XP_030844480.1.fas', 'XP_030828250.1.fas', 'XP_030840023.1.fas', 'XP_003728052.2.fas', 'XP_030836206.1.fas', 'XP_030838536.1.fas', 'NP_001295590.1.fas', 'XP_030830575.1.fas', 'XP_030856036.1.fas', 'XP_030833690.1.fas', 'XP_030845894.1.fas', 'XP_030839000.1.fas', 'XP_030827989.1.fas', 'XP_003725184.2.fas', 'XP_030848209.1.fas', 'XP_003725487.2.fas', 'XP_011664540.1.fas', 'XP_030833968.1.fas', 'XP_030852601.1.fas', 'XP_030838233.1.fas', 'XP_030828974.1.fas', 'XP_030841587.1.fas', 'XP_030855566.1.fas', 'XP_030830154.1.fas', 'XP_030850364.1.fas', 'XP_030832787.1.fas', 'XP_011677135.2.fas', 'XP_030828555.1.fas', 'XP_011674190.1.fas', 'XP_030835655.1.fas', 'XP_030840501.1.fas', 'NP_001154910.1.fas', 'XP_030837585.1.fas', 'XP_785947.2.fas', 'XP_030853615.1.fas', 'XP_030845524.1.fas', 'XP_030842929.1.fas', 'XP_030846156.1.fas', 'XP_030846455.1.fas', 'XP_795427.2.fas', 'XP_030830057.1.fas', 'XP_030836724.1.fas', 'XP_030830038.1.fas', 'XP_001189729.2.fas', 'XP_030841837.1.fas', 'XP_030831443.1.fas', 'XP_030842099.1.fas', 'XP_003728712.1.fas', 'NP_999759.1.fas', 'XP_030844214.1.fas', 'XP_799977.2.fas', 'XP_801771.1.fas', 'XP_030853982.1.fas', 'XP_784558.1.fas', 'XP_030831862.1.fas', 'XP_780332.4.fas', 'XP_790957.3.fas', 'XP_030853379.1.fas', 'XP_785736.3.fas', 'XP_030851557.1.fas', 'XP_030838959.1.fas', 'XP_030855469.1.fas', 'XP_787563.2.fas', 'XP_011662077.1.fas', 'XP_030841751.1.fas', 'XP_030843182.1.fas', 'XP_030832688.1.fas', 'XP_030834106.1.fas', 'NP_001229623.1.fas', 'XP_003724905.1.fas', 'XP_030849459.1.fas', 'XP_030847369.1.fas', 'XP_782821.1.fas', 'XP_794415.1.fas', 'XP_030831295.1.fas', 'XP_030840246.1.fas', 'XP_030850901.1.fas', 'XP_030855322.1.fas', 'XP_791119.2.fas', 'XP_788085.3.fas', 'XP_030830137.1.fas', 'XP_030854759.1.fas', 'XP_030840840.1.fas', 'XP_011677874.2.fas', 'XP_798605.2.fas', 'XP_011662982.2.fas', 'XP_003725530.1.fas', 'XP_003727918.1.fas', 'XP_780795.2.fas', 'XP_030837477.1.fas', 'XP_030846683.1.fas', 'XP_030829800.1.fas', 'XP_030847422.1.fas', 'XP_030837753.1.fas', 'XP_011670837.2.fas', 'XP_011673892.1.fas', 'XP_792619.4.fas', 'XP_030850993.1.fas', 'XP_030837250.1.fas', 'XP_030845341.1.fas', 'XP_030830015.1.fas', 'XP_784498.4.fas', 'XP_030832338.1.fas', 'XP_030838989.1.fas', 'XP_030856371.1.fas', 'XP_799247.3.fas', 'XP_030838372.1.fas', 'XP_003728808.3.fas', 'XP_030831820.1.fas', 'XP_030833408.1.fas', 'XP_003729701.2.fas', 'XP_030831725.1.fas', 'XP_030837974.1.fas', 'XP_030855894.1.fas', 'XP_030854614.1.fas', 'XP_030841273.1.fas', 'XP_798074.3.fas', 'XP_003726563.1.fas', 'XP_030851332.1.fas', 'XP_011661127.1.fas', 'XP_030846186.1.fas', 'XP_030846478.1.fas', 'XP_030855048.1.fas', 'XP_030831325.1.fas', 'XP_030837155.1.fas', 'XP_030848192.1.fas', 'XP_030854433.1.fas', 'NP_999731.2.fas', 'XP_030841557.1.fas', 'XP_030855692.1.fas', 'XP_030845242.1.fas', 'XP_030856156.1.fas', 'XP_003727279.2.fas', 'XP_030835430.1.fas', 'XP_030840764.1.fas', 'XP_030855303.1.fas', 'XP_030856272.1.fas', 'XP_011663849.2.fas', 'XP_800431.3.fas', 'XP_003725270.2.fas', 'XP_030852929.1.fas', 'XP_011661392.1.fas', 'XP_030840143.1.fas', 'XP_030832124.1.fas', 'XP_779985.5.fas', 'XP_011676613.2.fas', 'XP_030842255.1.fas', 'XP_030831371.1.fas', 'XP_030828608.1.fas', 'XP_030844526.1.fas', 'XP_030839332.1.fas', 'XP_011665926.2.fas', 'XP_030839587.1.fas', 'XP_011666911.1.fas', 'XP_030845216.1.fas', 'XP_011680778.1.fas', 'XP_030846961.1.fas', 'XP_796958.4.fas', 'XP_030833214.1.fas', 'XP_011671753.2.fas', 'XP_003730464.2.fas', 'XP_030845913.1.fas', 'XP_030850155.1.fas', 'XP_793296.1.fas', 'XP_030836115.1.fas', 'XP_030833433.1.fas', 'XP_030837090.1.fas', 'XP_030838920.1.fas', 'XP_030842356.1.fas', 'XP_030851541.1.fas', 'XP_030844425.1.fas', 'XP_011677695.2.fas', 'XP_030842950.1.fas', 'XP_786663.1.fas', 'XP_030851665.1.fas', 'XP_003729382.1.fas', 'XP_030852083.1.fas', 'XP_030833930.1.fas', 'XP_030843544.1.fas', 'XP_030837768.1.fas', 'XP_030839987.1.fas', 'XP_030846645.1.fas', 'XP_030848175.1.fas', 'XP_794275.4.fas', 'XP_011671039.2.fas', 'XP_030829854.1.fas', 'XP_003725891.2.fas', 'XP_030856295.1.fas', 'XP_784491.1.fas', 'XP_030838768.1.fas', 'XP_030831956.1.fas', 'XP_030842453.1.fas', 'XP_030841903.1.fas', 'XP_798717.1.fas', 'XP_030853497.1.fas', 'XP_030848838.1.fas', 'XP_030845082.1.fas', 'XP_030832006.1.fas', 'XP_030841421.1.fas', 'XP_030856396.1.fas', 'XP_030851560.1.fas', 'XP_030831253.1.fas', 'XP_030851941.1.fas', 'XP_030853097.1.fas', 'XP_030850839.1.fas', 'XP_030831474.1.fas', 'XP_011661152.1.fas', 'XP_030844223.1.fas', 'XP_785425.3.fas', 'XP_030855776.1.fas', 'XP_030841294.1.fas', 'XP_030845637.1.fas', 'XP_030850250.1.fas', 'XP_030836713.1.fas', 'XP_011662368.1.fas', 'XP_011671877.2.fas', 'XP_030839282.1.fas', 'XP_030828246.1.fas', 'XP_030830647.1.fas', 'XP_030856304.1.fas', 'XP_003723386.1.fas', 'XP_781936.3.fas', 'XP_030836134.1.fas', 'XP_003729076.1.fas', 'XP_001183758.2.fas', 'XP_030843203.1.fas', 'XP_030855431.1.fas', 'XP_030840056.1.fas', 'XP_030833172.1.fas', 'XP_030854690.1.fas', 'XP_793339.4.fas', 'XP_782356.2.fas', 'XP_003726778.2.fas', 'XP_794067.2.fas', 'XP_030841141.1.fas', 'XP_030847415.1.fas', 'XP_030837861.1.fas', 'XP_030845738.1.fas', 'XP_792912.2.fas', 'XP_030851227.1.fas', 'XP_011667242.1.fas', 'XP_011674581.2.fas', 'XP_030831935.1.fas', 'XP_011671257.1.fas', 'XP_030835790.1.fas', 'XP_030831630.1.fas', 'XP_792565.4.fas', 'XP_011668242.1.fas', 'XP_783565.2.fas', 'XP_030854896.1.fas', 'XP_030848731.1.fas', 'NP_999661.1.fas', 'XP_030835426.1.fas', 'XP_030855913.1.fas', 'XP_030830100.1.fas', 'XP_030841309.1.fas', 'XP_003725114.3.fas', 'XP_030834110.1.fas', 'XP_030843194.1.fas', 'NP_999829.1.fas', 'XP_011670922.2.fas', 'XP_011661039.2.fas', 'XP_011682852.2.fas', 'XP_030847613.1.fas', 'XP_783326.1.fas', 'XP_030830992.1.fas', 'XP_030846704.1.fas', 'XP_030831780.1.fas', 'XP_030828707.1.fas', 'XP_030854096.1.fas', 'XP_030829682.1.fas', 'XP_011669015.2.fas', 'XP_030853544.1.fas', 'XP_030855734.1.fas', 'XP_781189.4.fas', 'XP_030855237.1.fas', 'NP_001118232.1.fas', 'XP_030843106.1.fas', 'XP_011673610.1.fas', 'XP_030836088.1.fas', 'XP_030832128.1.fas', 'NP_001164287.1.fas', 'NP_001032400.1.fas', 'XP_030848137.1.fas', 'XP_003728827.2.fas', 'XP_030833553.1.fas', 'XP_789232.3.fas', 'XP_001181890.3.fas', 'XP_030829458.1.fas', 'XP_030849990.1.fas', 'XP_003731372.1.fas', 'XP_030835495.1.fas', 'XP_030837162.1.fas', 'XP_030847815.1.fas', 'XP_030842735.1.fas', 'XP_011681102.2.fas', 'XP_003726080.2.fas', 'XP_011661931.1.fas', 'XP_030854107.1.fas', 'XP_030841436.1.fas', 'XP_030852722.1.fas', 'XP_030845095.1.fas', 'XP_030838979.1.fas', 'XP_030847943.1.fas', 'XP_030837537.1.fas', 'XP_030842360.1.fas', 'XP_003725257.1.fas', 'XP_030849854.1.fas', 'XP_030851653.1.fas', 'XP_030833222.1.fas', 'XP_030836704.1.fas', 'NP_999813.1.fas', 'XP_003727115.1.fas', 'XP_011681236.1.fas', 'XP_030855146.1.fas', 'XP_030851518.1.fas', 'XP_783509.3.fas', 'XP_030836123.1.fas', 'XP_788119.4.fas', 'XP_001199227.2.fas', 'XP_030833622.1.fas', 'XP_030856210.1.fas', 'XP_783130.5.fas', 'XP_011682390.2.fas', 'XP_030833005.1.fas', 'XP_001187877.2.fas', 'XP_030828651.1.fas', 'XP_030848162.1.fas', 'XP_030830856.1.fas', 'XP_030855045.1.fas', 'XP_781306.4.fas', 'XP_030831644.1.fas', 'XP_030839607.1.fas', 'XP_030847745.1.fas', 'XP_795754.3.fas', 'XP_003729490.1.fas', 'XP_030847246.1.fas', 'XP_030846519.1.fas', 'XP_030840794.1.fas', 'XP_030837815.1.fas', 'XP_030836695.1.fas', 'XP_030851918.1.fas', 'NP_999832.1.fas', 'NP_001028182.1.fas', 'XP_800040.1.fas', 'XP_030827819.1.fas', 'XP_030836102.1.fas', 'XP_030853911.1.fas', 'XP_780963.2.fas', 'XP_030837087.1.fas', 'XP_030840003.1.fas', 'XP_030830989.1.fas', 'XP_030835654.1.fas', 'XP_030847962.1.fas', 'XP_030829961.1.fas', 'XP_030853686.1.fas', 'XP_030840091.1.fas', 'XP_030847364.1.fas', 'XP_030839820.1.fas', 'XP_030829699.1.fas', 'XP_011672270.1.fas', 'XP_030845996.1.fas', 'XP_030854657.1.fas', 'XP_030829540.1.fas', 'XP_030845693.1.fas', 'XP_030847764.1.fas', 'XP_786461.3.fas', 'XP_030831042.1.fas', 'XP_003725185.2.fas', 'XP_030832811.1.fas', 'XP_030844687.1.fas', 'XP_001183727.1.fas', 'XP_030836923.1.fas', 'XP_030840224.1.fas', 'XP_030830155.1.fas', 'XP_030840100.1.fas', 'XP_011678655.1.fas', 'NP_999818.2.fas', 'XP_011662487.1.fas', 'XP_003726598.2.fas', 'XP_030853871.1.fas', 'XP_011664397.1.fas', 'XP_030831991.1.fas', 'XP_011682122.1.fas', 'XP_003723523.2.fas', 'XP_786576.2.fas', 'XP_793412.2.fas', 'XP_030851878.1.fas', 'XP_030838175.1.fas', 'XP_030836940.1.fas', 'XP_780121.2.fas', 'XP_030842197.1.fas', 'XP_030844676.1.fas', 'XP_030838838.1.fas', 'XP_030845129.1.fas', 'XP_030830783.1.fas', 'XP_030841577.1.fas', 'XP_030828984.1.fas', 'XP_030829020.1.fas', 'XP_030854634.1.fas', 'XP_011667187.2.fas', 'XP_030850769.1.fas', 'XP_030849210.1.fas', 'XP_011660581.1.fas', 'XP_030837651.1.fas', 'XP_782249.1.fas', 'XP_030832053.1.fas', 'XP_030837575.1.fas', 'XP_030852898.1.fas', 'XP_011669763.1.fas', 'XP_791350.2.fas', 'XP_003724026.1.fas', 'XP_030846713.1.fas', 'XP_030829695.1.fas', 'XP_030856351.1.fas', 'XP_003725881.1.fas', 'XP_030851380.1.fas', 'XP_800714.4.fas', 'XP_030838588.1.fas', 'XP_011677259.1.fas', 'XP_030852128.1.fas', 'XP_030843111.1.fas', 'XP_030834526.1.fas', 'XP_030855693.1.fas', 'XP_030851733.1.fas', 'XP_030851230.1.fas', 'XP_030829820.1.fas', 'XP_011666204.2.fas', 'XP_030834968.1.fas', 'XP_011662336.2.fas', 'XP_030830233.1.fas', 'XP_030845167.1.fas', 'XP_030847749.1.fas', 'XP_030846934.1.fas', 'XP_030831092.1.fas', 'XP_030842448.1.fas', 'XP_030836664.1.fas', 'XP_030843033.1.fas', 'XP_796231.2.fas', 'XP_030837819.1.fas', 'XP_030836140.1.fas', 'XP_030833198.1.fas', 'XP_011660532.1.fas', 'XP_030830517.1.fas', 'XP_795397.1.fas', 'XP_030832339.1.fas', 'XP_030838554.1.fas', 'XP_030831492.1.fas', 'XP_030845946.1.fas', 'XP_030853071.1.fas', 'XP_030833241.1.fas', 'XP_030853572.1.fas', 'XP_030828116.1.fas', 'XP_030849837.1.fas', 'XP_791556.2.fas', 'XP_030837988.1.fas', 'XP_795931.3.fas', 'XP_030830959.1.fas', 'XP_030832571.1.fas', 'XP_030847625.1.fas', 'XP_030829626.1.fas', 'XP_787040.2.fas', 'XP_030843286.1.fas', 'XP_789273.2.fas', 'XP_030840015.1.fas', 'XP_003726383.1.fas', 'XP_011681201.1.fas', 'XP_030854009.1.fas', 'XP_030836417.1.fas', 'XP_787451.1.fas', 'XP_030846141.1.fas', 'XP_011661916.2.fas', 'XP_030844291.1.fas', 'XP_030850773.1.fas', 'XP_011675998.2.fas', 'XP_030838091.1.fas', 'XP_030837327.1.fas', 'XP_011660704.2.fas', 'XP_030854066.1.fas', 'XP_030842951.1.fas', 'XP_030830608.1.fas', 'NP_001121539.1.fas', 'XP_030832026.1.fas', 'XP_795786.2.fas', 'XP_030828209.1.fas', 'XP_030852810.1.fas', 'XP_030851862.1.fas', 'XP_030845278.1.fas', 'XP_011683169.2.fas', 'XP_011668560.2.fas', 'XP_798413.1.fas', 'XP_011662403.1.fas', 'XP_011663719.2.fas', 'XP_030840116.1.fas', 'XP_030853701.1.fas', 'XP_030828398.1.fas', 'XP_030844392.1.fas', 'XP_030855655.1.fas', 'XP_030851346.1.fas', 'XP_030829074.1.fas', 'XP_030852410.1.fas', 'XP_011668863.2.fas', 'XP_785837.1.fas', 'XP_030829956.1.fas', 'XP_001199084.2.fas', 'XP_011667766.2.fas', 'XP_030830562.1.fas', 'XP_011660547.1.fas', 'XP_030850251.1.fas', 'XP_784625.1.fas', 'XP_030852181.1.fas', 'XP_030841295.1.fas', 'XP_030851329.1.fas', 'XP_030843999.1.fas', 'XP_030855550.1.fas', 'XP_030840137.1.fas', 'XP_003725204.2.fas', 'XP_030853825.1.fas', 'XP_030835898.1.fas', 'XP_030832468.1.fas', 'XP_011677529.1.fas', 'XP_011665548.2.fas', 'XP_030836914.1.fas', 'XP_030831519.1.fas', 'XP_030840213.1.fas', 'XP_030850352.1.fas', 'XP_030834050.1.fas', 'XP_030831075.1.fas', 'XP_030832407.1.fas', 'XP_030854447.1.fas', 'XP_030844506.1.fas', 'XP_783052.4.fas', 'XP_011675165.2.fas', 'XP_030837441.1.fas', 'XP_030851401.1.fas', 'XP_030832464.1.fas', 'XP_030831934.1.fas', 'XP_791681.3.fas', 'NP_001354494.1.fas', 'XP_791382.3.fas', 'XP_030830749.1.fas', 'XP_783008.2.fas', 'XP_030851725.1.fas', 'XP_030831016.1.fas', 'XP_030846591.1.fas', 'XP_030853408.1.fas', 'XP_030855678.1.fas', 'XP_011682858.1.fas', 'XP_030845739.1.fas', 'XP_030832243.1.fas', 'XP_786407.1.fas', 'XP_030841308.1.fas', 'XP_030838142.1.fas', 'XP_030829186.1.fas', 'XP_011682732.1.fas', 'NP_999660.1.fas', 'XP_030856298.1.fas', 'XP_030830726.1.fas', 'XP_030834378.1.fas', 'XP_030840373.1.fas', 'XP_011680814.2.fas', 'XP_030850232.1.fas', 'XP_030830002.1.fas', 'XP_030851590.1.fas', 'XP_030853640.1.fas', 'XP_030833470.1.fas', 'XP_011673514.1.fas', 'XP_030830326.1.fas', 'XP_030838666.1.fas', 'XP_030845571.1.fas', 'XP_030830349.1.fas', 'XP_030841140.1.fas', 'XP_030854527.1.fas', 'XP_011672403.1.fas', 'XP_030855883.1.fas', 'XP_011665260.2.fas', 'XP_030846191.1.fas', 'XP_030828521.1.fas', 'XP_030853145.1.fas', 'XP_030833375.1.fas', 'XP_030850310.1.fas', 'XP_793338.3.fas', 'XP_030843004.1.fas', 'XP_030842258.1.fas', 'XP_784377.4.fas', 'XP_030836089.1.fas', 'XP_011671579.2.fas', 'XP_030837163.1.fas', 'XP_030832262.1.fas', 'XP_011675203.2.fas', 'XP_797311.2.fas', 'XP_030837563.1.fas', 'XP_011672121.1.fas', 'XP_011680486.2.fas', 'XP_030840019.1.fas', 'XP_782551.2.fas', 'XP_030852873.1.fas', 'XP_003725424.1.fas', 'XP_030831816.1.fas', 'XP_011680676.1.fas', 'XP_784339.3.fas', 'XP_030830395.1.fas', 'XP_030841843.1.fas', 'XP_030833219.1.fas', 'XP_011675603.2.fas', 'XP_030854622.1.fas', 'NP_999828.1.fas', 'XP_030847536.1.fas', 'XP_030832662.1.fas', 'XP_795548.3.fas', 'XP_011661670.2.fas', 'XP_030827794.1.fas', 'XP_030851668.1.fas', 'XP_003725061.2.fas', 'XP_030839118.1.fas', 'XP_030840352.1.fas', 'XP_030831879.1.fas', 'XP_030847680.1.fas', 'XP_011663387.2.fas', 'XP_030831282.1.fas', 'XP_011663981.2.fas', 'XP_030829986.1.fas', 'XP_030842409.1.fas', 'XP_781737.1.fas', 'XP_030832286.1.fas', 'XP_030852548.1.fas', 'XP_030845126.1.fas', 'XP_011681175.2.fas', 'XP_030845149.1.fas', 'XP_030841517.1.fas', 'XP_790671.4.fas', 'XP_030832214.1.fas', 'XP_030841330.1.fas', 'XP_011673506.2.fas', 'XP_030833391.1.fas', 'XP_030854852.1.fas', 'XP_802025.1.fas', 'XP_001197780.2.fas', 'XP_030842641.1.fas', 'NP_001091915.2.fas', 'XP_030847664.1.fas', 'XP_003727782.1.fas', 'XP_030827888.1.fas', 'XP_011674808.2.fas', 'XP_011666561.2.fas', 'XP_030851989.1.fas', 'XP_011675186.1.fas', 'XP_030829543.1.fas', 'XP_030854338.1.fas', 'XP_784312.1.fas', 'NP_999831.1.fas', 'XP_030828157.1.fas', 'XP_030855240.1.fas', 'XP_030843756.1.fas', 'XP_030856331.1.fas', 'XP_030838631.1.fas', 'XP_011662084.1.fas', 'XP_001178994.1.fas', 'XP_030835657.1.fas', 'XP_030855029.1.fas', 'XP_030838879.1.fas', 'XP_030839307.1.fas', 'XP_030833069.1.fas', 'XP_792047.2.fas', 'XP_030842866.1.fas', 'XP_003729190.1.fas', 'XP_792998.2.fas', 'XP_030828519.1.fas', 'XP_030830177.1.fas', 'XP_030842428.1.fas', 'XP_030835451.1.fas', 'XP_030828576.1.fas', 'XP_030835870.1.fas', 'XP_011666708.2.fas', 'XP_030835775.1.fas', 'XP_030848645.1.fas', 'XP_003729128.2.fas', 'XP_030836802.1.fas', 'XP_011664907.2.fas', 'XP_030855867.1.fas', 'XP_030832359.1.fas', 'XP_789887.2.fas', 'XP_030837058.1.fas', 'XP_011660676.1.fas', 'XP_030843274.1.fas', 'XP_030850463.1.fas', 'XP_784873.2.fas', 'XP_030839802.1.fas', 'XP_011677153.1.fas', 'XP_030830939.1.fas', 'XP_030855429.1.fas', 'XP_030847940.1.fas', 'XP_030829562.1.fas', 'XP_030833848.1.fas', 'XP_030827854.1.fas', 'XP_030839423.1.fas', 'XP_011667355.1.fas', 'XP_003730955.1.fas', 'XP_011671099.1.fas', 'XP_011684138.1.fas', 'XP_011675010.2.fas', 'XP_783506.3.fas', 'XP_030831727.1.fas', 'XP_030847923.1.fas', 'NP_001103718.1.fas', 'XP_011679692.1.fas', 'XP_030830630.1.fas', 'XP_791206.3.fas', 'XP_781347.3.fas', 'XP_030842392.1.fas', 'XP_011663594.2.fas', 'XP_030840366.1.fas', 'XP_030855701.1.fas', 'XP_030828416.1.fas', 'XP_030839608.1.fas', 'XP_030838875.1.fas', 'XP_001199744.2.fas', 'XP_030833642.1.fas', 'XP_796731.2.fas', 'XP_030844951.1.fas', 'XP_030850648.1.fas', 'XP_011680695.2.fas', 'XP_030827934.1.fas', 'XP_030831327.1.fas', 'XP_011660574.2.fas', 'XP_782165.2.fas', 'XP_003728880.2.fas', 'XP_030834904.1.fas', 'XP_030848190.1.fas', 'XP_011674798.1.fas', 'XP_030831624.1.fas', 'XP_030837454.1.fas', 'XP_030842569.1.fas', 'XP_782933.4.fas', 'XP_030828437.1.fas', 'XP_030833898.1.fas', 'XP_787834.3.fas', 'XP_003725284.1.fas', 'XP_799062.4.fas', 'XP_030845661.1.fas', 'XP_782335.4.fas', 'XP_030829520.1.fas', 'XP_030841753.1.fas', 'NP_001242990.1.fas', 'XP_030832972.1.fas', 'XP_030835582.1.fas', 'XP_011665852.2.fas', 'XP_783721.4.fas', 'XP_030827816.1.fas', 'XP_030854513.1.fas', 'XP_011682846.2.fas', 'XP_030828282.1.fas', 'XP_796062.4.fas', 'XP_785441.1.fas', 'XP_030856419.1.fas', 'XP_011673795.1.fas', 'XP_030836629.1.fas', 'XP_011660974.2.fas', 'XP_001197734.1.fas', 'XP_030841353.1.fas', 'XP_030854410.1.fas', 'NP_001268687.1.fas', 'XP_030832450.1.fas', 'XP_030831900.1.fas', 'XP_793565.3.fas', 'XP_011681487.2.fas', 'XP_030831992.1.fas', 'XP_011677583.1.fas', 'XP_003725684.1.fas', 'XP_030832218.1.fas', 'XP_001194806.2.fas', 'XP_784874.1.fas', 'XP_030847489.1.fas', 'XP_790682.2.fas', 'XP_030849060.1.fas', 'XP_030846628.1.fas', 'XP_791707.2.fas', 'XP_030851162.1.fas', 'XP_030831651.1.fas', 'XP_003728091.1.fas', 'XP_003724264.3.fas', 'XP_011667398.2.fas', 'XP_003727327.1.fas', 'XP_011665161.1.fas', 'XP_030845235.1.fas', 'XP_030830161.1.fas', 'XP_030836612.1.fas', 'XP_030843761.1.fas', 'XP_030835998.1.fas', 'XP_030832568.1.fas', 'XP_030853323.1.fas', 'XP_030836212.1.fas', 'XP_030843146.1.fas', 'XP_030838021.1.fas', 'XP_030837397.1.fas', 'XP_030831419.1.fas', 'XP_030835047.1.fas', 'XP_030846163.1.fas', 'XP_030853007.1.fas', 'XP_011680771.2.fas', 'XP_792693.3.fas', 'XP_030836886.1.fas', 'XP_011673131.2.fas', 'XP_030851943.1.fas', 'XP_030833482.1.fas', 'XP_030847350.1.fas', 'XP_030845433.1.fas', 'XP_030852384.1.fas', 'XP_030840115.1.fas', 'XP_792026.4.fas', 'XP_030838400.1.fas', 'XP_030832290.1.fas', 'XP_030840837.1.fas', 'XP_030832793.1.fas', 'XP_030835198.1.fas', 'XP_011682583.2.fas', 'XP_788478.2.fas', 'XP_030835409.1.fas', 'XP_030844905.1.fas', 'XP_011665322.2.fas', 'XP_011666679.1.fas', 'NP_001001768.2.fas', 'XP_030847156.1.fas', 'XP_011680518.2.fas', 'XP_030839330.1.fas', 'XP_797556.2.fas', 'XP_030828309.1.fas', 'NP_999767.2.fas', 'XP_030854465.1.fas', 'XP_780735.1.fas', 'XP_011670345.2.fas', 'XP_030831457.1.fas', 'XP_011673413.2.fas', 'XP_030834255.1.fas', 'XP_030842657.1.fas', 'XP_030845530.1.fas', 'XP_030838627.1.fas', 'XP_030828141.1.fas', 'XP_011676459.2.fas', 'XP_030845317.1.fas', 'XP_030830043.1.fas', 'XP_786841.2.fas', 'XP_030838560.1.fas', 'XP_011663384.2.fas', 'XP_030831887.1.fas', 'XP_784774.4.fas', 'XP_030855111.1.fas', 'XP_030840075.1.fas', 'XP_030852870.1.fas', 'XP_030831213.1.fas', 'XP_030843849.1.fas', 'NP_001229581.1.fas', 'XP_011662895.2.fas', 'XP_030834236.1.fas', 'XP_030841162.1.fas', 'XP_030850879.1.fas', 'XP_003728798.1.fas', 'XP_030827800.1.fas', 'XP_030835097.1.fas', 'XP_030852152.1.fas', 'XP_030832661.1.fas', 'XP_030850282.1.fas', 'XP_001177619.1.fas', 'XP_030850781.1.fas', 'XP_030856228.1.fas', 'XP_030830268.1.fas', 'XP_001185164.2.fas', 'XP_030852973.1.fas', 'XP_011680078.2.fas', 'XP_791745.2.fas', 'XP_030846369.1.fas', 'XP_030841646.1.fas', 'XP_030854221.1.fas', 'XP_030843095.1.fas', 'XP_030837747.1.fas', 'XP_003724778.1.fas', 'XP_779912.1.fas', 'XP_030854827.1.fas', 'XP_011666485.2.fas', 'XP_030833376.1.fas', 'XP_030855336.1.fas', 'XP_001199154.2.fas', 'XP_030830801.1.fas', 'XP_030855012.1.fas', 'XP_030833052.1.fas', 'XP_030856247.1.fas', 'XP_030836439.1.fas', 'XP_030839877.1.fas', 'XP_030852754.1.fas', 'XP_030832599.1.fas', 'XP_011661519.2.fas', 'XP_030832564.1.fas', 'XP_030843868.1.fas', 'XP_030832945.1.fas', 'XP_030831116.1.fas', 'XP_030853099.1.fas', 'XP_030833254.1.fas', 'XP_030842379.1.fas', 'XP_030828724.1.fas', 'XP_030832008.1.fas', 'XP_030836155.1.fas', 'XP_030830325.1.fas', 'XP_030855911.1.fas', 'XP_030850934.1.fas', 'NP_001189420.1.fas', 'XP_003725140.2.fas', 'XP_030849592.1.fas', 'XP_030834386.1.fas', 'XP_030834980.1.fas', 'XP_030848114.1.fas', 'XP_030842779.1.fas', 'XP_800733.1.fas', 'XP_003727126.2.fas', 'XP_030835792.1.fas', 'XP_030843996.1.fas', 'XP_030850084.1.fas', 'XP_011676653.2.fas', 'XP_030839056.1.fas', 'XP_030855378.1.fas', 'XP_011664218.2.fas', 'XP_030855983.1.fas', 'XP_030838085.1.fas', 'XP_785533.4.fas', 'XP_011675377.2.fas', 'XP_030831440.1.fas', 'XP_030854571.1.fas', 'XP_030847366.1.fas', 'XP_030834242.1.fas', 'XP_030844430.1.fas', 'XP_785120.1.fas', 'XP_011684052.2.fas', 'XP_011677612.2.fas', 'XP_030830976.1.fas', 'XP_030829609.1.fas', 'XP_791245.3.fas', 'XP_030836403.1.fas', 'XP_011683636.2.fas', 'XP_030847309.1.fas', 'XP_003728783.1.fas', 'XP_001181022.2.fas', 'XP_030836500.1.fas', 'XP_030828371.1.fas', 'XP_030834928.1.fas', 'XP_030833525.1.fas', 'XP_030837485.1.fas', 'NP_999659.1.fas', 'XP_790476.3.fas', 'XP_001179750.2.fas', 'XP_784072.2.fas', 'XP_030849671.1.fas', 'XP_030841632.1.fas', 'XP_030853787.1.fas', 'XP_030852301.1.fas', 'XP_030854171.1.fas', 'XP_030829266.1.fas', 'XP_780940.2.fas', 'XP_030841158.1.fas', 'XP_011675816.1.fas', 'XP_003728973.2.fas', 'XP_011670235.1.fas', 'XP_030835389.1.fas', 'XP_030839769.1.fas', 'XP_030837986.1.fas', 'XP_011661325.2.fas', 'XP_030830075.1.fas', 'XP_030850746.1.fas', 'XP_030842029.1.fas', 'XP_030846174.1.fas', 'XP_030855763.1.fas', 'XP_785580.4.fas', 'XP_030851651.1.fas', 'XP_030831162.1.fas', 'XP_030853581.1.fas', 'XP_030846418.1.fas', 'XP_030851352.1.fas', 'XP_030854377.1.fas', 'XP_030854689.1.fas', 'XP_030838986.1.fas', 'XP_030852223.1.fas', 'XP_001198416.2.fas', 'XP_799248.3.fas', 'XP_030844933.1.fas', 'XP_030852004.1.fas', 'XP_030842446.1.fas', 'XP_030849077.1.fas', 'XP_030852623.1.fas', 'XP_011665533.2.fas', 'XP_030851175.1.fas', 'XP_030835218.1.fas', 'XP_003731093.1.fas', 'XP_030833095.1.fas', 'XP_011672089.1.fas', 'XP_030833504.1.fas', 'XP_780699.5.fas', 'XP_030832182.1.fas', 'NP_001075434.1.fas', 'XP_030836521.1.fas', 'XP_030840704.1.fas', 'XP_030827939.1.fas', 'XP_030844659.1.fas', 'XP_792216.3.fas', 'XP_011665488.1.fas', 'XP_030851838.1.fas', 'XP_030830475.1.fas', 'XP_011669206.2.fas', 'XP_030831593.1.fas', 'XP_030836963.1.fas', 'NP_999689.1.fas', 'XP_795259.3.fas', 'XP_003724405.1.fas', 'XP_030833099.1.fas', 'XP_030855548.1.fas', 'XP_011678174.2.fas', 'XP_030829727.1.fas', 'XP_030852640.1.fas', 'XP_030851731.1.fas', 'XP_030839041.1.fas', 'XP_011674143.1.fas', 'XP_030847400.1.fas', 'XP_003724267.2.fas', 'XP_030835983.1.fas', 'XP_030848092.1.fas', 'XP_030837556.1.fas', 'XP_030839765.1.fas', 'XP_030853338.1.fas', 'XP_030855148.1.fas', 'XP_011678779.1.fas', 'XP_030834124.1.fas', 'XP_011675933.2.fas', 'XP_798456.4.fas', 'XP_011662113.2.fas', 'XP_030831193.1.fas', 'XP_781184.2.fas', 'XP_030851584.1.fas', 'XP_030834869.1.fas', 'XP_030834894.1.fas', 'XP_030845565.1.fas', 'XP_030830282.1.fas', 'XP_030855569.1.fas', 'XP_789669.4.fas', 'XP_030835759.1.fas', 'XP_030830879.1.fas', 'XP_030843082.1.fas', 'XP_030832276.1.fas', 'XP_030852545.1.fas', 'XP_011662750.1.fas', 'XP_011669248.2.fas', 'XP_030829422.1.fas', 'XP_003724897.1.fas', 'XP_030842107.1.fas', 'XP_779905.1.fas', 'XP_003729504.2.fas', 'XP_030832870.1.fas', 'XP_003730511.2.fas', 'NP_999655.1.fas', 'XP_030840746.1.fas', 'XP_030835412.1.fas', 'XP_030849978.1.fas', 'XP_030856174.1.fas', 'XP_781558.2.fas', 'XP_786955.2.fas', 'XP_011666361.1.fas', 'XP_030833045.1.fas', 'XP_030851158.1.fas', 'XP_003723829.1.fas', 'XP_030840462.1.fas', 'XP_030830713.1.fas', 'XP_030835833.1.fas', 'XP_030854659.1.fas', 'XP_790579.3.fas', 'XP_030855721.1.fas', 'XP_003725957.2.fas', 'XP_030848522.1.fas', 'XP_030840062.1.fas', 'XP_011675854.1.fas', 'XP_030835635.1.fas', 'XP_780304.2.fas', 'XP_790710.4.fas', 'XP_011676863.2.fas', 'XP_030839744.1.fas', 'XP_030841188.1.fas', 'XP_030841175.1.fas', 'XP_030829306.1.fas', 'XP_787941.2.fas', 'XP_030846159.1.fas', 'XP_030834105.1.fas', 'XP_030832676.1.fas', 'XP_030838518.1.fas', 'XP_030855276.1.fas', 'XP_011679732.2.fas', 'XP_030850253.1.fas', 'XP_003726289.2.fas', 'XP_030841297.1.fas', 'NP_001232902.1.fas', 'XP_030840036.1.fas', 'XP_030835661.1.fas', 'XP_030853322.1.fas', 'XP_011662242.2.fas', 'XP_030844495.1.fas', 'XP_011670825.1.fas', 'XP_030838304.1.fas', 'XP_030845510.1.fas', 'XP_003727726.1.fas', 'XP_030839815.1.fas', 'XP_030847957.1.fas', 'XP_030845358.1.fas', 'XP_030840986.1.fas', 'XP_030831174.1.fas', 'XP_030842774.1.fas', 'XP_030831955.1.fas', 'XP_003724612.2.fas', 'XP_011664353.1.fas', 'XP_011661485.2.fas', 'XP_030854445.1.fas', 'XP_030831589.1.fas', 'XP_030853197.1.fas', 'XP_030844620.1.fas', 'XP_030831574.1.fas', 'XP_011679150.1.fas', 'XP_030834052.1.fas', 'XP_030847475.1.fas', 'XP_030850653.1.fas', 'XP_784821.1.fas', 'XP_030850956.1.fas', 'XP_030855676.1.fas', 'XP_030848751.1.fas', 'XP_030854129.1.fas', 'XP_030845110.1.fas', 'XP_030836537.1.fas', 'XP_782480.2.fas', 'XP_798395.4.fas', 'XP_030835762.1.fas', 'XP_030840085.1.fas', 'XP_003724835.2.fas', 'XP_030853692.1.fas', 'XP_792448.4.fas', 'XP_030837001.1.fas', 'XP_030841403.1.fas', 'XP_785826.3.fas', 'XP_780706.3.fas', 'XP_030833285.1.fas', 'XP_030842071.1.fas', 'XP_003726534.1.fas', 'XP_030834170.1.fas', 'XP_030851963.1.fas', 'XP_011676134.2.fas', 'XP_030849764.1.fas', 'XP_030840935.1.fas', 'XP_003730697.1.fas', 'XP_030845910.1.fas', 'XP_003731211.1.fas', 'XP_030831439.1.fas', 'XP_030837093.1.fas', 'XP_030847919.1.fas', 'XP_030848557.1.fas', 'XP_011682758.2.fas', 'XP_011669285.1.fas', 'XP_030833430.1.fas', 'XP_030837590.1.fas', 'XP_030848273.1.fas', 'XP_030833314.1.fas', 'XP_783905.2.fas', 'XP_030846364.1.fas', 'XP_030828699.1.fas', 'XP_784737.3.fas', 'NP_001123279.1.fas', 'XP_030831974.1.fas', 'XP_030852317.1.fas', 'XP_001190869.1.fas', 'XP_030842755.1.fas', 'XP_783071.4.fas', 'XP_030841327.1.fas', 'XP_030838190.1.fas', 'XP_030854845.1.fas', 'NP_001091928.1.fas', 'XP_011670549.1.fas', 'XP_030841921.1.fas', 'XP_030831555.1.fas', 'XP_030831536.1.fas', 'XP_030841942.1.fas', 'XP_001200412.3.fas', 'XP_786231.3.fas', 'XP_030832763.1.fas', 'XP_030837462.1.fas', 'XP_030846368.1.fas', 'XP_011669689.1.fas', 'XP_030854104.1.fas', 'XP_030832144.1.fas', 'NP_001268690.1.fas', 'XP_030850535.1.fas', 'XP_030836588.1.fas', 'XP_030853865.1.fas', 'XP_030844529.1.fas', 'XP_030829887.1.fas', 'XP_030835720.1.fas', 'XP_030833550.1.fas', 'XP_030849368.1.fas', 'XP_030845870.1.fas', 'XP_030836954.1.fas', 'NP_001091924.1.fas', 'XP_030828523.1.fas', 'XP_030854596.1.fas', 'XP_782232.1.fas', 'XP_030839838.1.fas', 'XP_030854095.1.fas', 'XP_030855234.1.fas', 'XP_030828123.1.fas', 'XP_030850211.1.fas', 'XP_030830522.1.fas', 'XP_011668544.1.fas', 'XP_030840939.1.fas', 'XP_780142.3.fas', 'XP_030835568.1.fas', 'XP_030837645.1.fas', 'XP_030831435.1.fas', 'XP_001198442.2.fas', 'XP_782050.2.fas', 'XP_011676083.1.fas', 'XP_030829916.1.fas', 'XP_030852871.1.fas', 'XP_030828295.1.fas', 'XP_030828325.1.fas', 'XP_030843905.1.fas', 'XP_030855032.1.fas', 'XP_780230.2.fas', 'XP_030830227.1.fas', 'XP_030848231.1.fas', 'XP_011679652.2.fas', 'XP_030853465.1.fas', 'XP_030847484.1.fas', 'XP_030843524.1.fas', 'XP_030845754.1.fas', 'XP_792865.3.fas', 'XP_030829487.1.fas', 'XP_011674337.2.fas', 'XP_787993.3.fas', 'XP_011668901.2.fas', 'XP_030834532.1.fas', 'XP_030844946.1.fas', 'XP_030841666.1.fas', 'XP_030839554.1.fas', 'XP_030842130.1.fas', 'XP_030851224.1.fas', 'XP_792838.1.fas', 'XP_030850085.1.fas', 'XP_003723774.1.fas', 'XP_030832198.1.fas', 'XP_030834315.1.fas', 'XP_794287.4.fas', 'XP_030831330.1.fas', 'XP_001200128.2.fas', 'XP_030828190.1.fas', 'XP_792019.1.fas', 'XP_030828493.1.fas', 'XP_011667944.1.fas', 'XP_030847934.1.fas', 'XP_030831835.1.fas', 'NP_001116993.1.fas', 'XP_030841142.1.fas', 'XP_030850417.1.fas', 'XP_030830324.1.fas', 'XP_030830627.1.fas', 'XP_003728603.2.fas', 'XP_003726133.2.fas', 'XP_784173.4.fas', 'XP_030845952.1.fas', 'XP_784941.1.fas', 'XP_030847426.1.fas', 'XP_030827913.1.fas', 'XP_030835487.1.fas', 'XP_030842403.1.fas', 'XP_030851214.1.fas', 'XP_011681413.2.fas', 'XP_030855090.1.fas', 'XP_030847807.1.fas', 'XP_030842727.1.fas', 'XP_030842224.1.fas', 'NP_001020380.1.fas', 'XP_030828315.1.fas', 'XP_030853272.1.fas', 'XP_030830217.1.fas', 'XP_030854187.1.fas', 'XP_030840741.1.fas', 'XP_030833960.1.fas', 'XP_030850303.1.fas', 'XP_030830133.1.fas', 'XP_030853977.1.fas', 'XP_030831897.1.fas', 'XP_030836199.1.fas', 'XP_030845543.1.fas', 'XP_011660937.1.fas', 'XP_001199278.3.fas', 'XP_030836240.1.fas', 'XP_011671722.2.fas', 'XP_030848302.1.fas', 'XP_011661190.1.fas', 'XP_030836829.1.fas', 'XP_030839467.1.fas', 'XP_030850894.1.fas', 'XP_030854934.1.fas', 'XP_030851614.1.fas', 'XP_030832372.1.fas', 'XP_030839743.1.fas', 'XP_030851530.1.fas', 'XP_030846881.1.fas', 'NP_001229591.1.fas', 'XP_030842297.1.fas', 'XP_030828637.1.fas', 'XP_011663545.1.fas', 'XP_795058.4.fas', 'XP_030855901.1.fas', 'XP_030829195.1.fas', 'XP_030855604.1.fas', 'XP_030848723.1.fas', 'XP_785475.4.fas', 'XP_792408.1.fas', 'XP_030835782.1.fas', 'XP_785252.4.fas', 'XP_030842205.1.fas', 'XP_030856219.1.fas', 'XP_795313.4.fas', 'XP_030836808.1.fas', 'XP_030827831.1.fas', 'XP_783327.3.fas', 'XP_003726567.1.fas', 'XP_800559.3.fas', 'XP_011672918.2.fas', 'XP_030842306.1.fas', 'XP_011670752.1.fas', 'XP_030850406.1.fas', 'XP_030845562.1.fas', 'XP_797451.4.fas', 'XP_030836867.1.fas', 'XP_030853074.1.fas', 'XP_003730434.2.fas', 'XP_030853577.1.fas', 'XP_011663592.2.fas', 'XP_030850722.1.fas', 'XP_003728555.2.fas', 'XP_030845326.1.fas', 'XP_030851639.1.fas', 'XP_030845920.1.fas', 'XP_030828170.1.fas', 'XP_030830655.1.fas', 'XP_791983.2.fas', 'XP_030847946.1.fas', 'XP_030837532.1.fas', 'XP_030828738.1.fas', 'XP_030847643.1.fas', 'XP_030830339.1.fas', 'XP_003725854.1.fas', 'XP_030844732.1.fas', 'XP_003727182.1.fas', 'XP_030845349.1.fas', 'XP_030842860.1.fas', 'XP_030832117.1.fas', 'XP_030840448.1.fas', 'XP_030839602.1.fas', 'XP_787702.2.fas', 'XP_030841614.1.fas', 'XP_003730833.1.fas', 'XP_030836907.1.fas', 'XP_030848740.1.fas', 'XP_030843556.1.fas', 'XP_003727013.1.fas', 'XP_030845225.1.fas', 'XP_030840124.1.fas', 'XP_030828654.1.fas', 'XP_030848167.1.fas', 'XP_030830255.1.fas', 'NP_001075433.1.fas', 'XP_030838658.1.fas', 'XP_011668383.1.fas', 'XP_011677174.1.fas', 'XP_784296.4.fas', 'XP_030836195.1.fas', 'XP_030852706.1.fas', 'XP_030840069.1.fas', 'XP_001201921.2.fas', 'XP_797758.4.fas', 'XP_011674208.2.fas', 'XP_030832612.1.fas', 'XP_030845696.1.fas', 'XP_030833206.1.fas', 'XP_030835588.1.fas', 'XP_030832680.1.fas', 'XP_030855461.1.fas', 'XP_030840006.1.fas', 'XP_030830689.1.fas', 'XP_790505.2.fas', 'XP_030836623.1.fas', 'NP_999631.1.fas', 'XP_030835476.1.fas', 'XP_030853436.1.fas', 'XP_030828970.1.fas', 'XP_797510.4.fas', 'XP_030856234.1.fas', 'XP_030838831.1.fas', 'XP_030846375.1.fas', 'XP_030848146.1.fas', 'XP_030842228.1.fas', 'XP_030842744.1.fas', 'XP_030831660.1.fas', 'XP_011679842.1.fas', 'XP_030850528.1.fas', 'XP_030847761.1.fas', 'XP_783059.2.fas', 'XP_030835838.1.fas', 'XP_030837237.1.fas', 'XP_003726730.3.fas', 'XP_784157.5.fas', 'XP_011668159.1.fas', 'XP_030841336.1.fas', 'XP_030844915.1.fas', 'NP_999665.1.fas', 'XP_030831582.1.fas', 'XP_780532.5.fas', 'XP_030854294.1.fas', 'XP_030851496.1.fas', 'XP_030842281.1.fas', 'XP_030833075.1.fas', 'XP_030850010.1.fas', 'XP_030839972.1.fas', 'XP_030853729.1.fas', 'XP_797512.2.fas', 'XP_030835769.1.fas', 'XP_030852575.1.fas', 'XP_030841362.1.fas', 'XP_030830468.1.fas', 'XP_001199035.1.fas', 'XP_011670202.2.fas', 'XP_791482.4.fas', 'XP_030850658.1.fas', 'XP_003723457.1.fas', 'XP_784522.2.fas', 'XP_030840719.1.fas', 'XP_796543.1.fas', 'XP_011672621.1.fas', 'XP_030831013.1.fas', 'XP_011683104.1.fas', 'XP_030831737.1.fas', 'XP_030837547.1.fas', 'XP_030837044.1.fas', 'XP_030847636.1.fas', 'XP_011661817.1.fas', 'XP_030854021.1.fas', 'XP_030829336.1.fas', 'XP_030841145.1.fas', 'XP_030828197.1.fas', 'XP_011671279.2.fas', 'XP_030834135.1.fas', 'XP_780272.1.fas', 'XP_003725432.3.fas', 'XP_030832646.1.fas', 'XP_030847512.1.fas', 'XP_030830007.1.fas', 'XP_786763.2.fas', 'XP_797673.1.fas', 'XP_030835900.1.fas', 'XP_030837465.1.fas', 'XP_030830868.1.fas', 'XP_030837166.1.fas', 'XP_030828391.1.fas', 'XP_030830293.1.fas', 'XP_030849906.1.fas', 'XP_030844665.1.fas', 'XP_030833919.1.fas', 'XP_011665347.2.fas', 'XP_011681721.2.fas', 'XP_785185.2.fas', 'XP_011680157.1.fas', 'XP_030850315.1.fas', 'XP_030850616.1.fas', 'XP_030854248.1.fas', 'XP_030830201.1.fas', 'XP_784671.4.fas', 'XP_781951.4.fas', 'XP_030830702.1.fas', 'XP_030846603.1.fas', 'XP_011662304.2.fas', 'XP_030833770.1.fas', 'XP_030845974.1.fas', 'XP_030836850.1.fas', 'XP_030833273.1.fas', 'XP_030842579.1.fas', 'XP_030841829.1.fas', 'XP_003725099.2.fas', 'XP_011679079.1.fas', 'XP_003726231.2.fas', 'XP_030838065.1.fas', 'XP_790473.2.fas', 'XP_030830026.1.fas', 'XP_030855835.1.fas', 'XP_011662704.2.fas', 'XP_791913.4.fas', 'XP_030846700.1.fas', 'XP_003726479.2.fas', 'XP_011683328.2.fas', 'XP_011676777.2.fas', 'XP_030837065.1.fas', 'NP_999728.1.fas', 'XP_030840338.1.fas', 'XP_030833282.1.fas', 'XP_030839412.1.fas', 'XP_030842588.1.fas', 'XP_030831451.1.fas', 'XP_030828128.1.fas', 'XP_030849260.1.fas', 'XP_030842651.1.fas', 'XP_011682730.2.fas', 'XP_030842352.1.fas', 'XP_030829618.1.fas', 'XP_782443.3.fas', 'XP_030828898.1.fas', 'XP_788933.2.fas', 'XP_030851029.1.fas', 'XP_003724182.2.fas', 'XP_030843161.1.fas', 'XP_001176098.2.fas', 'NP_999724.1.fas', 'XP_030846144.1.fas', 'XP_780968.3.fas', 'XP_030837494.1.fas', 'XP_030855589.1.fas', 'XP_792401.3.fas', 'XP_030835744.1.fas', 'XP_030856222.1.fas', 'XP_030841595.1.fas', 'XP_785092.3.fas', 'XP_011681023.1.fas', 'XP_030830445.1.fas', 'XP_030832795.1.fas', 'XP_789594.2.fas', 'XP_030838406.1.fas', 'XP_030843062.1.fas', 'XP_782308.3.fas', 'XP_030856194.1.fas', 'XP_030849965.1.fas', 'XP_030832802.1.fas', 'XP_030847777.1.fas', 'XP_030837406.1.fas', 'XP_030844021.1.fas', 'XP_030853626.1.fas', 'XP_030855456.1.fas', 'NP_999705.1.fas', 'XP_030846466.1.fas', 'XP_030830064.1.fas', 'XP_030843140.1.fas', 'XP_030844724.1.fas', 'XP_011677594.2.fas', 'XP_780135.1.fas', 'XP_011666377.2.fas', 'XP_011684062.2.fas', 'XP_781846.2.fas', 'XP_030839717.1.fas', 'XP_030828509.1.fas', 'XP_030827947.1.fas', 'XP_030829256.1.fas', 'XP_030847255.1.fas', 'XP_030836530.1.fas', 'XP_030829239.1.fas', 'XP_030832490.1.fas', 'XP_030844092.1.fas', 'XP_030828341.1.fas', 'XP_030837448.1.fas', 'XP_792420.3.fas', 'XP_780096.3.fas', 'XP_801707.2.fas', 'XP_030835441.1.fas', 'XP_030846066.1.fas', 'XP_011680377.2.fas', 'XP_030855727.1.fas', 'XP_030839409.1.fas', 'XP_030835517.1.fas', 'XP_030846433.1.fas', 'XP_799857.1.fas', 'XP_030838072.1.fas', 'XP_030829392.1.fas', 'XP_003725951.2.fas', 'XP_030854514.1.fas', 'XP_030847303.1.fas', 'XP_011673840.2.fas', 'XP_030851531.1.fas', 'XP_011667183.2.fas', 'XP_030832670.1.fas', 'XP_030847027.1.fas', 'XP_785762.1.fas', 'XP_011676649.1.fas', 'XP_030836828.1.fas', 'XP_030837472.1.fas', 'XP_785956.3.fas', 'XP_030829700.1.fas', 'XP_030837853.1.fas', 'XP_030832270.1.fas', 'XP_003725380.1.fas', 'XP_030833961.1.fas', 'XP_030845765.1.fas', 'XP_797656.4.fas', 'XP_030842193.1.fas', 'XP_030848200.1.fas', 'XP_011677579.1.fas', 'XP_030854485.1.fas', 'XP_030836598.1.fas', 'XP_030833043.1.fas', 'XP_030842749.1.fas', 'XP_030838853.1.fas', 'XP_011673267.2.fas', 'XP_011666303.2.fas', 'XP_798540.1.fas', 'XP_030841152.1.fas', 'XP_030835985.1.fas', 'XP_030854036.1.fas', 'XP_030832575.1.fas', 'XP_030851337.1.fas', 'XP_030853519.1.fas', 'XP_011663002.2.fas', 'XP_791068.4.fas', 'XP_030843134.1.fas', 'XP_030829097.1.fas', 'XP_030847094.1.fas', 'XP_781182.2.fas', 'XP_030829594.1.fas', 'XP_782915.4.fas', 'XP_011665010.1.fas', 'NP_999658.2.fas', 'XP_030828512.1.fas', 'XP_011679947.2.fas', 'XP_030833940.1.fas', 'XP_793816.2.fas', 'XP_003729242.1.fas', 'XP_030842795.1.fas', 'XP_030836544.1.fas', 'XP_030832419.1.fas', 'XP_030843313.1.fas', 'XP_030832476.1.fas', 'XP_030841552.1.fas', 'XP_786731.1.fas', 'XP_782458.4.fas', 'XP_011679509.2.fas', 'XP_011682666.2.fas', 'XP_030845228.1.fas', 'XP_030832857.1.fas', 'XP_030844653.1.fas', 'XP_030837289.1.fas', 'XP_030855609.1.fas', 'XP_011661042.1.fas', 'XP_030846639.1.fas', 'XP_011661698.1.fas', 'XP_030842764.1.fas', 'XP_030834698.1.fas', 'XP_030842861.1.fas', 'XP_003726043.2.fas', 'XP_030839994.1.fas', 'XP_030848844.1.fas', 'XP_003726367.2.fas', 'NP_001029121.1.fas', 'XP_030841379.1.fas', 'XP_011674149.1.fas', 'XP_030839548.1.fas', 'XP_030835456.1.fas', 'XP_003724563.1.fas', 'XP_030838912.1.fas', 'XP_030828853.1.fas', 'XP_030839897.1.fas', 'XP_030855798.1.fas', 'XP_001196473.3.fas', 'XP_030841287.1.fas', 'XP_030829565.1.fas', 'XP_030836792.1.fas', 'XP_030854977.1.fas', 'XP_030829342.1.fas', 'XP_780929.2.fas', 'XP_030849172.1.fas', 'XP_030836148.1.fas', 'XP_030828377.1.fas', 'XP_030833020.1.fas', 'XP_030830275.1.fas', 'XP_030853816.1.fas', 'XP_030855344.1.fas', 'XP_030840723.1.fas', 'XP_030832887.1.fas', 'XP_030832281.1.fas', 'XP_003724568.2.fas', 'XP_030849374.1.fas', 'XP_030847147.1.fas', 'XP_003731160.2.fas', 'XP_030839506.1.fas', 'XP_030839927.1.fas', 'XP_030830719.1.fas', 'XP_030832434.1.fas', 'XP_030848128.1.fas', 'XP_030828318.1.fas', 'XP_030831661.1.fas', 'XP_030840495.1.fas', 'XP_030839321.1.fas', 'XP_030845992.1.fas', 'XP_003729462.2.fas', 'XP_798530.3.fas', 'XP_030845369.1.fas', 'XP_030838580.1.fas', 'XP_030839721.1.fas', 'XP_030829998.1.fas', 'XP_030847966.1.fas', 'XP_011681071.2.fas', 'XP_780716.3.fas', 'XP_011677983.2.fas', 'XP_030842943.1.fas', 'XP_030838636.1.fas', 'XP_030837580.1.fas', 'XP_781220.4.fas', 'XP_788322.2.fas', 'XP_030828274.1.fas', 'XP_011676940.1.fas', 'XP_030838011.1.fas', 'XP_030831429.1.fas', 'XP_030832979.1.fas', 'NP_999733.1.fas', 'XP_030836824.1.fas', 'XP_011669947.2.fas', 'XP_783728.2.fas', 'XP_011680279.2.fas', 'XP_030854523.1.fas', 'XP_030850180.1.fas', 'XP_784207.2.fas', 'XP_030843892.1.fas', 'XP_030833474.1.fas', 'XP_781040.1.fas', 'XP_795649.3.fas', 'XP_030845076.1.fas', 'XP_030836870.1.fas', 'XP_799798.1.fas', 'XP_030838546.1.fas', 'XP_030840971.1.fas', 'XP_030851948.1.fas', 'XP_030846926.1.fas', 'XP_003724415.1.fas', 'XP_030846623.1.fas', 'XP_030850636.1.fas', 'XP_030852018.1.fas', 'XP_030847482.1.fas', 'XP_790097.4.fas', 'NP_999699.1.fas', 'XP_030853463.1.fas', 'XP_030830494.1.fas', 'XP_030847410.1.fas', 'XP_030854704.1.fas', 'XP_030829110.1.fas', 'XP_030855984.1.fas', 'XP_030830848.1.fas', 'XP_030831336.1.fas', 'XP_030834915.1.fas', 'XP_003724179.1.fas', 'XP_030829737.1.fas', 'XP_011670594.2.fas', 'XP_801793.1.fas', 'XP_011670569.2.fas', 'XP_011681798.1.fas', 'XP_011675990.2.fas', 'XP_030844708.1.fas', 'XP_011675695.2.fas', 'XP_030835501.1.fas', 'XP_003726551.1.fas', 'XP_011668223.2.fas', 'NP_999729.1.fas', 'XP_030850285.1.fas', 'XP_030836238.1.fas', 'XP_030848981.1.fas', 'XP_011673250.2.fas', 'XP_788825.3.fas', 'XP_030855184.1.fas', 'XP_030851024.1.fas', 'XP_030828790.1.fas', 'XP_030852271.1.fas', 'XP_030834515.1.fas', 'XP_791368.1.fas', 'XP_783407.1.fas', 'XP_030854820.1.fas', 'XP_030830791.1.fas', 'XP_030841565.1.fas', 'XP_030833539.1.fas', 'XP_783517.4.fas', 'XP_030831614.1.fas', 'XP_011674256.1.fas', 'XP_790458.1.fas', 'XP_787312.1.fas', 'XP_030835823.1.fas', 'XP_797599.4.fas', 'XP_787202.4.fas', 'XP_011672390.1.fas', 'XP_011668949.2.fas', 'XP_030830147.1.fas', 'XP_030845710.1.fas', 'XP_030850377.1.fas', 'XP_030830444.1.fas', 'XP_030840830.1.fas', 'XP_030838104.1.fas', 'XP_011672400.2.fas', 'XP_030840735.1.fas', 'XP_786004.2.fas', 'XP_030853206.1.fas', 'XP_030855076.1.fas', 'XP_030855575.1.fas', 'XP_011676913.2.fas', 'XP_030853294.1.fas', 'XP_001185403.2.fas', 'XP_780603.3.fas', 'XP_030832121.1.fas', 'XP_011674038.2.fas', 'XP_030842856.1.fas', 'XP_030845158.1.fas', 'XP_011670198.2.fas', 'XP_001199778.2.fas', 'XP_011663358.1.fas', 'XP_030839013.1.fas', 'XP_030837723.1.fas', 'XP_030838196.1.fas', 'XP_030829152.1.fas', 'XP_030854561.1.fas', 'XP_011676216.2.fas', 'XP_030846808.1.fas', 'XP_030836182.1.fas', 'XP_030829388.1.fas', 'XP_030832521.1.fas', 'XP_030842955.1.fas', 'XP_011664579.2.fas', 'XP_030840335.1.fas', 'XP_030855251.1.fas', 'XP_030842018.1.fas', 'XP_785190.3.fas', 'XP_011678393.2.fas', 'XP_030837596.1.fas', 'XP_030855476.1.fas', 'XP_787272.1.fas', 'XP_030854443.1.fas', 'XP_030847757.1.fas', 'XP_030845179.1.fas', 'XP_030856290.1.fas', 'XP_793925.1.fas', 'XP_030837807.1.fas', 'XP_030854264.1.fas', 'XP_011668964.2.fas', 'XP_030850339.1.fas', 'XP_030841300.1.fas', 'XP_030840786.1.fas', 'XP_030842155.1.fas', 'XP_030828508.1.fas', 'XP_011680376.2.fas', 'XP_794360.1.fas', 'XP_030836615.1.fas', 'XP_791550.3.fas', 'XP_011667999.2.fas', 'XP_030855554.1.fas', 'XP_030839687.1.fas', 'XP_030839379.1.fas', 'XP_030830741.1.fas', 'XP_795937.2.fas', 'XP_003730497.1.fas', 'XP_030852185.1.fas', 'XP_030837996.1.fas', 'XP_030845632.1.fas', 'XP_030830098.1.fas', 'XP_030848654.1.fas', 'XP_030836813.1.fas', 'XP_030842039.1.fas', 'XP_030845015.1.fas', 'XP_011679113.2.fas', 'XP_780969.4.fas', 'XP_011668429.2.fas', 'XP_030833114.1.fas', 'XP_003724812.2.fas', 'XP_030842372.1.fas', 'XP_030834888.1.fas', 'XP_030832003.1.fas', 'XP_003728609.2.fas', 'XP_781241.2.fas', 'XP_030852233.1.fas', 'XP_786170.4.fas', 'XP_030841805.1.fas', 'XP_030837302.1.fas', 'XP_003726010.1.fas', 'XP_030850739.1.fas', 'XP_030852414.1.fas', 'XP_791229.1.fas', 'XP_785753.4.fas', 'XP_003724546.2.fas', 'XP_011665372.2.fas', 'XP_030849335.1.fas', 'XP_030832751.1.fas', 'XP_030844955.1.fas', 'XP_030839547.1.fas', 'XP_030842123.1.fas', 'XP_030855397.1.fas', 'XP_030836966.1.fas', 'XP_011661846.2.fas', 'XP_030841319.1.fas', 'XP_030850320.1.fas', 'XP_030844018.1.fas', 'XP_011674000.2.fas', 'XP_030835817.1.fas', 'XP_030836547.1.fas', 'XP_030829069.1.fas', 'XP_011670679.2.fas', 'XP_030839128.1.fas', 'XP_003725552.2.fas', 'XP_785018.4.fas', 'XP_030830337.1.fas', 'XP_030843213.1.fas', 'XP_030847948.1.fas', 'XP_011666392.2.fas', 'XP_030853352.1.fas', 'XP_030835683.1.fas', 'XP_003728171.1.fas', 'XP_030831826.1.fas', 'XP_011681533.2.fas', 'XP_030850168.1.fas', 'XP_011662184.2.fas', 'XP_003726296.2.fas', 'XP_030847506.1.fas', 'XP_030830432.1.fas', 'XP_030855324.1.fas', 'XP_030833399.1.fas', 'XP_030828530.1.fas', 'XP_030827983.1.fas', 'XP_030851279.1.fas', 'XP_030836566.1.fas', 'XP_030856255.1.fas', 'XP_003730333.2.fas', 'NP_001020382.1.fas', 'XP_030840199.1.fas', 'XP_011680506.2.fas', 'XP_030834921.1.fas', 'XP_030831302.1.fas', 'XP_030854117.1.fas', 'XP_030835478.1.fas', 'XP_779900.1.fas', 'XP_003729002.2.fas', 'XP_796285.2.fas', 'XP_030847127.1.fas', 'XP_030854233.1.fas', 'XP_030832557.1.fas', 'XP_003725435.1.fas', 'XP_003723645.1.fas', 'XP_793277.1.fas', 'XP_030836109.1.fas', 'XP_030830379.1.fas', 'XP_030839844.1.fas', 'XP_030830384.1.fas', 'XP_030852767.1.fas', 'XP_030830982.1.fas', 'XP_030847906.1.fas', 'XP_030829524.1.fas', 'XP_030839465.1.fas', 'XP_787944.2.fas', 'XP_030832976.1.fas', 'XP_030851387.1.fas', 'XP_030838572.1.fas', 'XP_011681770.1.fas', 'XP_030831790.1.fas', 'XP_788678.4.fas', 'XP_011682026.1.fas', 'XP_030853916.1.fas', 'XP_030841181.1.fas', 'XP_030856335.1.fas', 'XP_030843857.1.fas', 'XP_030844159.1.fas', 'XP_030842329.1.fas', 'XP_791272.1.fas', 'XP_030843175.1.fas', 'XP_030829528.1.fas', 'NP_999730.1.fas', 'XP_003724366.1.fas', 'NP_001001474.1.fas', 'XP_030844212.1.fas', 'XP_030839105.1.fas', 'XP_030850808.1.fas', 'XP_030845694.1.fas', 'XP_003724104.2.fas', 'XP_030854353.1.fas', 'XP_030847544.1.fas', 'XP_011673820.2.fas', 'XP_784294.4.fas', 'XP_030841334.1.fas', 'XP_030851873.1.fas', 'XP_030834563.1.fas', 'XP_030838183.1.fas', 'XP_030827972.1.fas', 'XP_011672389.1.fas', 'XP_030849971.1.fas', 'XP_030835118.1.fas', 'XP_030853887.1.fas', 'XP_030829865.1.fas', 'XP_030839322.1.fas', 'XP_011663006.1.fas', 'XP_030847809.1.fas', 'XP_789660.4.fas', 'XP_030836006.1.fas', 'XP_030839094.1.fas', 'NP_999633.1.fas', 'XP_030853137.1.fas', 'XP_011682899.1.fas', 'XP_030835489.1.fas', 'XP_030830451.1.fas', 'XP_030844985.1.fas', 'XP_030841658.1.fas', 'XP_030834142.1.fas', 'XP_030836791.1.fas', 'XP_030835539.1.fas', 'XP_030839427.1.fas', 'XP_785435.3.fas', 'XP_011672256.1.fas', 'XP_030856386.1.fas', 'XP_011661096.2.fas', 'XP_030851073.1.fas', 'XP_030829947.1.fas', 'XP_003727830.1.fas', 'XP_030839292.1.fas', 'XP_030836124.1.fas', 'XP_030830354.1.fas', 'XP_003723396.1.fas', 'XP_030856314.1.fas', 'XP_030836427.1.fas', 'XP_030853015.1.fas', 'XP_030851338.1.fas', 'XP_030850846.1.fas', 'XP_030838530.1.fas', 'XP_011667674.2.fas', 'XP_030830257.1.fas', 'XP_011682499.1.fas', 'XP_003724390.2.fas', 'XP_030848165.1.fas', 'XP_030830851.1.fas', 'XP_030850343.1.fas', 'XP_030841387.1.fas', 'XP_030854877.1.fas', 'NP_999680.1.fas', 'XP_030840793.1.fas', 'XP_030837717.1.fas', 'XP_780900.1.fas', 'XP_030847466.1.fas', 'XP_030837812.1.fas', 'NP_001032719.1.fas', 'XP_030836692.1.fas', 'XP_030851448.1.fas', 'XP_030830805.1.fas', 'XP_030829882.1.fas', 'XP_030847287.1.fas', 'XP_030836073.1.fas', 'XP_030845875.1.fas', 'NP_001165523.1.fas', 'XP_003729276.1.fas', 'XP_030838199.1.fas', 'XP_030855934.1.fas', 'XP_030830127.1.fas', 'XP_785479.2.fas', 'XP_030854823.1.fas', 'XP_030849026.1.fas', 'XP_030832141.1.fas', 'XP_030855084.1.fas', 'XP_003726775.1.fas', 'XP_030844040.1.fas', 'XP_030851905.1.fas', 'XP_030854625.1.fas', 'XP_030855858.1.fas', 'XP_003729730.2.fas', 'XP_011666337.2.fas', 'XP_011671183.2.fas', 'XP_011678546.2.fas', 'XP_030838640.1.fas', 'XP_030831883.1.fas', 'XP_030835923.1.fas', 'XP_030831285.1.fas', 'XP_030855416.1.fas', 'XP_030853365.1.fas', 'XP_030855115.1.fas', 'XP_030850214.1.fas', 'XP_001185328.2.fas', 'XP_030833772.1.fas', 'XP_030835502.1.fas', 'XP_030842085.1.fas', 'XP_782803.3.fas', 'XP_030850812.1.fas', 'XP_030828568.1.fas', 'XP_030844646.1.fas', 'XP_003731137.2.fas', 'XP_030851827.1.fas', 'XP_003725625.1.fas', 'XP_030830497.1.fas', 'XP_030851105.1.fas', 'XP_030852956.1.fas', 'XP_003727421.2.fas', 'XP_030843697.1.fas', 'XP_030841044.1.fas', 'XP_030830721.1.fas', 'XP_782003.3.fas', 'XP_030853744.1.fas', 'XP_030838867.1.fas', 'XP_787104.4.fas', 'XP_030828004.1.fas', 'XP_030830927.1.fas', 'XP_011679985.1.fas', 'XP_030838661.1.fas', 'XP_030828825.1.fas', 'XP_030828404.1.fas', 'XP_789627.1.fas', 'XP_003727499.1.fas', 'XP_030832347.1.fas', 'XP_003724453.2.fas', 'XP_030830594.1.fas', 'XP_788426.3.fas', 'XP_030849826.1.fas', 'XP_030845938.1.fas', 'XP_798099.3.fas', 'XP_793697.2.fas', 'XP_794646.2.fas', 'XP_030831830.1.fas', 'XP_791493.1.fas', 'XP_030844461.1.fas', 'XP_030855458.1.fas', 'XP_030829932.1.fas', 'XP_011679974.1.fas', 'XP_030853348.1.fas', 'XP_793436.3.fas', 'XP_030829654.1.fas', 'XP_030831171.1.fas', 'XP_030844225.1.fas', 'XP_786851.4.fas', 'XP_030839431.1.fas', 'XP_003729456.2.fas', 'XP_783545.1.fas', 'XP_011670360.2.fas', 'XP_030840911.1.fas', 'XP_030854608.1.fas', 'XP_800280.2.fas', 'XP_030828199.1.fas', 'XP_030827829.1.fas', 'XP_788839.1.fas', 'XP_030855770.1.fas', 'XP_030844749.1.fas', 'XP_011667761.2.fas', 'XP_030838602.1.fas', 'XP_030839882.1.fas', 'XP_030835367.1.fas', 'XP_802079.1.fas', 'XP_030855157.1.fas', 'XP_011670001.1.fas', 'XP_789451.3.fas', 'XP_001201542.3.fas', 'XP_030832248.1.fas', 'XP_030843542.1.fas', 'XP_782052.1.fas', 'XP_030841391.1.fas', 'XP_030839684.1.fas', 'XP_030829738.1.fas', 'XP_003729413.1.fas', 'XP_030832492.1.fas', 'XP_030835862.1.fas', 'XP_030845484.1.fas', 'XP_001195486.3.fas', 'XP_790898.3.fas', 'XP_030832103.1.fas', 'XP_781983.4.fas', 'XP_030845187.1.fas', 'XP_784879.3.fas', 'XP_003724912.2.fas', 'XP_030844002.1.fas', 'XP_030832227.1.fas', 'XP_030851844.1.fas', 'XP_030847470.1.fas', 'XP_030842455.1.fas', 'XP_030855252.1.fas', 'XP_030842518.1.fas', 'XP_003724180.2.fas', 'XP_030843460.1.fas', 'XP_003729615.1.fas', 'XP_030841406.1.fas', 'XP_030847676.1.fas', 'XP_030828808.1.fas', 'XP_030847973.1.fas', 'XP_030837507.1.fas', 'XP_030853369.1.fas', 'XP_030829052.1.fas', 'XP_030850218.1.fas', 'XP_011678302.2.fas', 'XP_030829551.1.fas', 'XP_030856095.1.fas', 'XP_011683763.1.fas', 'NP_001138965.1.fas', 'XP_030831453.1.fas', 'XP_011663889.2.fas', 'XP_030839637.1.fas', 'XP_030831377.1.fas', 'XP_003729248.3.fas', 'XP_786298.1.fas', 'XP_786766.1.fas', 'XP_780015.1.fas', 'XP_001189437.1.fas', 'XP_030851263.1.fas', 'XP_011666173.2.fas', 'XP_030854745.1.fas', 'XP_030845282.1.fas', 'XP_001199075.1.fas', 'XP_787234.2.fas', 'XP_030855351.1.fas', 'XP_030839082.1.fas', 'XP_030832892.1.fas', 'XP_030844395.1.fas', 'XP_003728244.2.fas', 'XP_030843344.1.fas', 'XP_030852380.1.fas', 'XP_030835245.1.fas', 'XP_001198150.3.fas', 'XP_794433.1.fas', 'XP_030853650.1.fas', 'XP_030842397.1.fas', 'XP_030836445.1.fas', 'XP_030838676.1.fas', 'XP_001193210.3.fas', 'XP_030853589.1.fas', 'XP_030854681.1.fas', 'XP_003726205.2.fas', 'XP_030846181.1.fas', 'XP_030843280.1.fas', 'XP_030841150.1.fas', 'XP_030842903.1.fas', 'XP_030832577.1.fas', 'XP_030852244.1.fas', 'XP_782185.3.fas', 'XP_781236.2.fas', 'XP_030842421.1.fas', 'XP_790028.2.fas', 'XP_030847825.1.fas', 'XP_003728283.2.fas', 'XP_030831924.1.fas', 'XP_030852347.1.fas', 'XP_030839908.1.fas', 'XP_030834993.1.fas', 'XP_030844019.1.fas', 'XP_030853753.1.fas', 'XP_030830111.1.fas', 'XP_030856151.1.fas', 'XP_011664299.2.fas', 'XP_030848223.1.fas', 'XP_030833209.1.fas', 'XP_030851912.1.fas', 'XP_030837070.1.fas', 'XP_030831703.1.fas', 'XP_030856338.1.fas', 'XP_030849431.1.fas', 'XP_030842922.1.fas', 'XP_791458.1.fas', 'XP_030831869.1.fas', 'XP_030838657.1.fas', 'XP_030836167.1.fas', 'XP_011664894.1.fas', 'XP_030855102.1.fas', 'XP_030829996.1.fas', 'XP_030828215.1.fas', 'XP_030834193.1.fas', 'XP_030851980.1.fas', 'XP_011680223.2.fas', 'XP_001179358.2.fas', 'XP_030835515.1.fas', 'XP_030853271.1.fas', 'XP_011681371.1.fas', 'XP_030847896.1.fas', 'XP_030828615.1.fas', 'XP_030842248.1.fas', 'XP_780700.2.fas', 'XP_784245.3.fas', 'XP_030845443.1.fas', 'XP_030836099.1.fas', 'XP_030852309.1.fas', 'XP_001198703.3.fas', 'NP_999651.1.fas', 'XP_011678514.1.fas', 'XP_030830130.1.fas', 'XP_011675887.2.fas', 'XP_782504.1.fas', 'XP_030841655.1.fas', 'XP_003727230.2.fas', 'XP_779901.1.fas', 'XP_030827910.1.fas', 'XP_030839946.1.fas', 'XP_030846981.1.fas', 'XP_030848149.1.fas', 'XP_030842724.1.fas', 'XP_011675154.2.fas', 'XP_030834943.1.fas', 'XP_030836095.1.fas', 'XP_011677089.1.fas', 'XP_003724322.2.fas', 'XP_786658.2.fas', 'XP_030851274.1.fas', 'XP_030856181.1.fas', 'XP_001193519.3.fas', 'XP_030832211.1.fas', 'XP_030854752.1.fas', 'XP_030841335.1.fas', 'XP_011676549.2.fas', 'XP_030855940.1.fas', 'XP_030830153.1.fas', 'XP_011682093.2.fas', 'XP_030854838.1.fas', 'XP_030848261.1.fas', 'XP_030856237.1.fas', 'XP_030832459.1.fas', 'XP_796151.2.fas', 'XP_011666022.1.fas', 'XP_030845123.1.fas', 'XP_011669164.2.fas', 'XP_030835252.1.fas', 'XP_003726841.2.fas', 'XP_030846151.1.fas', 'XP_030835576.1.fas', 'XP_030832986.1.fas', 'XP_030830050.1.fas', 'NP_999834.1.fas', 'XP_001180189.1.fas', 'XP_001195605.3.fas', 'XP_030843189.1.fas', 'XP_030833422.1.fas', 'XP_030844158.1.fas', 'XP_030842328.1.fas', 'XP_794681.2.fas', 'XP_030853985.1.fas', 'XP_011673227.2.fas', 'XP_011672959.2.fas', 'XP_030837510.1.fas', 'XP_030842644.1.fas', 'XP_030856082.1.fas', 'XP_030833794.1.fas', 'XP_011666599.2.fas', 'NP_001001475.1.fas', 'XP_799859.1.fas', 'XP_030831147.1.fas', 'XP_003727900.3.fas', 'XP_030835157.1.fas', 'XP_030851739.1.fas', 'XP_030844958.1.fas', 'XP_030840805.1.fas', 'XP_784123.3.fas', 'XP_030833921.1.fas', 'XP_030853233.1.fas', 'XP_011682695.2.fas', 'XP_030833500.1.fas', 'XP_011660580.2.fas', 'XP_030828354.1.fas', 'XP_030842709.1.fas', 'XP_030850598.1.fas', 'XP_797532.4.fas', 'XP_030836026.1.fas', 'XP_030831947.1.fas', 'XP_030852922.1.fas', 'XP_792159.3.fas', 'XP_003727156.2.fas', 'XP_779971.3.fas', 'XP_011672159.2.fas', 'XP_030832230.1.fas', 'XP_783375.1.fas', 'XP_030829946.1.fas', 'XP_030845093.1.fas', 'XP_790472.4.fas', 'XP_011668776.2.fas', 'XP_030839125.1.fas', 'XP_030839426.1.fas', 'XP_796223.3.fas', 'XP_030849254.1.fas', 'XP_030833822.1.fas', 'XP_787968.2.fas', 'XP_030839868.1.fas', 'XP_781203.4.fas', 'XP_030840527.1.fas', 'XP_030838910.1.fas', 'XP_003723464.2.fas', 'XP_030854500.1.fas', 'NP_001229579.1.fas', 'XP_030835648.1.fas', 'XP_780864.3.fas', 'XP_030847911.1.fas', 'XP_030842631.1.fas', 'XP_782054.2.fas', 'XP_030830549.1.fas', 'XP_030829030.1.fas', 'XP_030835000.1.fas', 'XP_030846124.1.fas', 'XP_011669438.1.fas', 'XP_797966.3.fas', 'XP_030854348.1.fas', 'XP_030852138.1.fas', 'XP_003726232.2.fas', 'XP_030830526.1.fas', 'XP_784350.3.fas', 'XP_030845556.1.fas', 'XP_011681597.2.fas', 'XP_030836171.1.fas', 'XP_030830126.1.fas', 'XP_011675891.2.fas', 'XP_011678001.1.fas', 'XP_030842484.1.fas', 'XP_030840754.1.fas', 'XP_030836571.1.fas', 'XP_784446.1.fas', 'XP_030846600.1.fas', 'XP_030853799.1.fas', 'XP_030855514.1.fas', 'XP_030837165.1.fas', 'XP_030846692.1.fas', 'XP_030850982.1.fas', 'XP_030851807.1.fas', 'XP_030854727.1.fas', 'XP_030840973.1.fas', 'XP_030855712.1.fas', 'NP_999798.1.fas', 'XP_030840375.1.fas', 'XP_030836872.1.fas', 'XP_011665623.1.fas', 'XP_001197240.2.fas', 'XP_030835397.1.fas', 'XP_030837047.1.fas', 'XP_030832062.1.fas', 'XP_030839872.1.fas', 'XP_030852252.1.fas', 'XP_030842915.1.fas', 'XP_030828169.1.fas', 'XP_030831113.1.fas', 'XP_030830595.1.fas', 'XP_030837998.1.fas', 'XP_030850582.1.fas', 'XP_001198564.2.fas', 'XP_030832161.1.fas', 'XP_030852652.1.fas', 'XP_030841045.1.fas', 'XP_011678498.2.fas', 'XP_030835892.1.fas', 'XP_030839689.1.fas', 'XP_011667462.1.fas', 'XP_030844942.1.fas', 'XP_030832843.1.fas', 'XP_030851723.1.fas', 'XP_796639.3.fas', 'XP_030855312.1.fas', 'XP_030854794.1.fas', 'XP_030836689.1.fas', 'XP_030830825.1.fas', 'XP_003730305.2.fas', 'XP_030838866.1.fas', 'XP_030829259.1.fas', 'NP_001027542.1.fas', 'XP_030830223.1.fas', 'XP_787044.3.fas', 'XP_030837148.1.fas', 'XP_030853726.1.fas', 'XP_030854209.1.fas', 'XP_030856124.1.fas', 'XP_011675025.1.fas', 'XP_030853402.1.fas', 'XP_003730849.1.fas', 'XP_030837792.1.fas', 'XP_030832820.1.fas', 'XP_011677197.2.fas', 'XP_030832725.1.fas', 'XP_797071.2.fas', 'XP_030834974.1.fas', 'XP_011676930.2.fas', 'XP_030854441.1.fas', 'XP_780035.1.fas', 'XP_030834372.1.fas', 'XP_030849065.1.fas', 'XP_030847572.1.fas', 'XP_798916.1.fas', 'XP_011663576.2.fas', 'XP_030854666.1.fas', 'XP_030832325.1.fas', 'XP_030831757.1.fas', 'XP_790077.1.fas', 'XP_030835366.1.fas', 'XP_786427.2.fas', 'XP_030829339.1.fas', 'XP_030830640.1.fas', 'XP_030836430.1.fas', 'XP_003726457.2.fas', 'XP_030850470.1.fas', 'XP_011666116.1.fas', 'XP_030837392.1.fas', 'XP_011678443.1.fas', 'XP_030837969.1.fas', 'XP_030855874.1.fas', 'XP_030845333.1.fas', 'XP_030855889.1.fas', 'XP_030854247.1.fas', 'XP_030838497.1.fas', 'XP_030838169.1.fas', 'XP_011679175.1.fas', 'XP_030846528.1.fas', 'XP_030842176.1.fas', 'XP_011668341.2.fas', 'XP_792075.2.fas', 'XP_030853890.1.fas', 'XP_011670464.2.fas', 'XP_030847774.1.fas', 'XP_011668893.1.fas', 'XP_011676911.2.fas', 'XP_030839335.1.fas', 'XP_030843943.1.fas', 'XP_030831319.1.fas', 'XP_030837497.1.fas', 'XP_011675926.1.fas', 'XP_030845135.1.fas', 'XP_030849909.1.fas', 'XP_786322.2.fas', 'XP_030850388.1.fas', 'XP_030836411.1.fas', 'XP_030829318.1.fas', 'XP_001180150.3.fas', 'XP_030841495.1.fas', 'XP_030842938.1.fas', 'XP_030838622.1.fas', 'XP_011682551.1.fas', 'XP_030830362.1.fas', 'XP_030849128.1.fas', 'XP_030851529.1.fas', 'XP_797937.1.fas', 'XP_030840013.1.fas', 'XP_030838506.1.fas', 'XP_030844891.1.fas', 'XP_011668020.1.fas', 'XP_003723289.2.fas', 'XP_030842519.1.fas', 'XP_030835560.1.fas', 'XP_030831452.1.fas', 'XP_030847553.1.fas', 'XP_030852134.1.fas', 'XP_030854344.1.fas', 'XP_030832304.1.fas', 'XP_030852816.1.fas', 'XP_030837005.1.fas', 'XP_030832020.1.fas', 'XP_011677462.1.fas']

	fasta_file_list = fasta_files[0:859]
	#fasta_file_list = fasta_files[859:1718]
	#fasta_file_list = fasta_files[1718:2579]
	Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_file_list)

	#clean_up_iqtree_files()
	#subset_boot_file()
	#edit_tree_files("loci.treefile","single_locus_trees.nwk")
	#edit_tree_files("loci.boottrees", "single_locus_trees_boot_subset.nwk")
	#clean_gene_trees("single_locus_trees_boot_subset.nwk", "clean_trees.nwk")
	#root_trees()

if __name__ == "__main__":
	main()
