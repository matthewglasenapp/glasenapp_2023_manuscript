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
vcf_file = "/hb/scratch/mglasena/data/genotypes/strongylocentrotidae/insertions_removed.vcf.gz"

# Bed file containing a record for each protein coding gene in the S. purpuratus assembly. See the ncbi/ directory for scripts to generate this file
protein_coding_genes_bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_genes/protein_coding_genes.bed"

# Directory containing the output files from mrna.py
bed_file_dir = "/hb/scratch/mglasena/mrna_cov/"

# Sample bed file from the output of mrna.py
bed_file = "/hb/scratch/mglasena/mrna_cov/pallidus_SRR5767285.regions.bed.gz"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
# Strongylocentrotidae Subset
subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Franciscanus outgroup
#subset_sample_list = ['droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']

# Specify thresholds for filtering. 
min_cov_threshold = 10

prop_1x_threshold = 0.75

prop_10x_threshold = 0.75

# Required gap between two adjacent loci in base pairs 
required_gap = 20000

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

# Populate subset_mean_coverage_spur5_exons variable with the samples included in the subset_sample_list variable
def subset_coverage_dict():
	for sample in subset_sample_list:
		subset_mean_coverage_spur5_exons[sample] = mean_coverage_spur5_exons[sample]

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
	for i in range(4):
		for sample in subset_sample_list:
			header.append(sample)
	
	writer.writerow(header)

	# Write row for each mRNA containing mean depth, prop1x, prop10x, and prop20x for each species/sample
	counter = 0
	for key,value in rna_dict.items():
		row = [key] + value[0] + value[1] + value[2] + value[3]
		writer.writerow(row)
		counter += 1

	csv_file.close()
	
	print("There were {} mRNA records pre-filter".format(len(rna_dict)))
	print("{} mRNA records written to all_rna.csv".format(counter))

# Filter rna_dict by coverage metrics. Add mRNA's passing filter to passed_rna_dict
def filter_rna_dict():
	record_counter = 0
	failed_counter = 0
	passed_counter = 0

	for key, value in rna_dict.items():
		record_counter += 1
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
				passed_counter += 1 
			else:
				failed_counter += 1
		else:
			failed_counter += 1

	print("{} mRNA records processed".format(record_counter))
	print("{} mRNA records passed intial coverage depth filters".format(passed_counter))
	print("{} mRNA records failed intial coverage depth filters".format(failed_counter))
	print("{} mRNA records in passed_rna_dict".format(len(passed_rna_dict)))

# Get file of mRNA records from gff3 file
def get_mrna_gff():
	get_mrna_from_gff = '''awk '$3 == "mRNA"' {} > mrna_records.txt'''.format(gff_file)
	os.system(get_mrna_from_gff)
	
# Create dictionary for mRNAs passing initial filter in the format of "Scaffold":[[scaffold,rna,start,stop]]. This dictionary will be used for filtering by distance on chromosome
def create_scaffold_dict():
	
	# Creates mapping of rna: parent_gene for mRNAs in passed_rna_dict 
	mt_rna_counter = 0
	record_counter = 0

	with open("mrna_records.txt", "r") as f:
		for line in f:
			rna = line.split("\t")[8].split(";")[0].split("rna-")[1]
			if rna in passed_rna_dict.keys():
				scaffold = line.split("\t")[0]
				start = line.split("\t")[3]
				stop = line.split("\t")[4]
				parent_gene = line.split("\t")[8].split(";")[1].split("gene-")[1]
				
				if scaffold == "NC_001453.1":
					mt_rna_counter += 1

				else:
					mrna_gene_dict[rna] = parent_gene

					if scaffold_dict.get(scaffold):
						scaffold_dict[scaffold].append([scaffold,rna,start,stop])
					else:
						scaffold_dict[scaffold] = [[scaffold,rna,start,stop]]

					record_counter += 1

	os.system("rm mrna_records.txt")
	print("{} mitochondrial mRNAs were removed".format(mt_rna_counter))
	print("{} mRNA records in scaffold_dict".format(record_counter))
	print("{} mRNA records in mrna_gene_dict".format(len(mrna_gene_dict)))

# Filter mRNA records in scaffold_dict so that all remaining mRNAs are >= required_gap apart from each other
def check_proximity():
	filter_counter = 0
	
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
				rna_list.pop(rna_counter)
				filter_counter +=1 
			
			else:
				last_stop = int(stop)
				rna_counter +=1
	
	record_counter = 0
	for value in scaffold_dict.values():
		for item in value:
			record_counter += 1

	print("{} mRNAs removed to satisfy gap filter".format(filter_counter))
	print("{} mRNAs records remain in scaffold_dict".format(record_counter))
	print("{} mRNAs passed both coverage depth and gap filters".format(record_counter))

# Get list of passed mRNAs and create new filtered_mrna_gene_dict that only includes mRNAs that passed the gap filter. 
def get_passed_rnas():
	for rna_list in scaffold_dict.values():
		for rna in rna_list:
			passed_rnas.append(rna[1])

	for key,value in mrna_gene_dict.items():
		if key in passed_rnas:
			filtered_mrna_gene_dict[key] = value

	print("{} mRNA records in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))

# Write new bed file ("unlinked_loci.bed") of parent genes with an mRNA transcript that passed all filters. This will be used to build the initial vcf2fasta alignments 
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

# Vcf2fasta makes alignments for all isoforms of each gene. This function identifies the unwanted isoform alignment files and deletes them.
def remove_redundant_isoforms():
	# List of redundant isoforms to delete following vcf2fasta
	records_to_delete = []

	get_cds_gff = '''awk '$3 == "CDS"' sco_gff.gff > cds.gff'''
	os.system(get_cds_gff)

	records = open("cds.gff","r").read().splitlines()

	for record in records:
		cds_name = record.split("\t")[8].split(";")[0].split("cds-")[1]
		parent_rna_name = record.split("\t")[8].split(";")[1].split("rna-")[1]
		cds_parent_rna_dict[cds_name] = parent_rna_name
	
	passed_rnas_lst = list(filtered_mrna_gene_dict.keys())

	for key,value in cds_parent_rna_dict.items():
		if not value in passed_rnas_lst:
			records_to_delete.append(key)

	for record in records_to_delete:
		cds_parent_rna_dict.pop(record)
		delete = "rm vcf2fasta_CDS/{}.fas".format(record)
		os.system(delete)

	#os.system("rm sco_gff.gff")
	#os.system("rm cds.gff")

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_CDS/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_CDS/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

# Run iqturee on the fasta files of mRNAs that passed all filters
def identify_no_variant_no_parsimony():
	run_iqtree = "iqtree2 -S vcf2fasta_CDS/ -m TESTONLY --prefix loci -T 8"
	os.system(run_iqtree)

def remove_no_variant_no_parsimony():
	os.mkdir('no_variant_no_parsimony')
	no_variant_no_parsimony_lst = []
	
	get_no_variant_no_parsimony = '''cat loci.log | grep -E "No parsimony|No variant" > no_variant_no_parsimony.txt'''

	os.system(get_no_variant_no_parsimony)

	with open("no_variant_no_parsimony.txt", "r") as f:
		for line in f:
			no_variant_no_parsimony_lst.append(line.split(" ")[6].strip())
	
	for file in no_variant_no_parsimony_lst:
		cds = file.split(".fas")[0]
		rna = cds_parent_rna_dict[cds]
		passed_rnas.remove(rna)
		filtered_mrna_gene_dict.pop(rna)
		cds_parent_rna_dict.pop(cds)

		
		move = "mv vcf2fasta_CDS/{} no_variant_no_parsimony/".format(file)
		os.system(move)

	os.system("rm no_variant_no_parsimony.txt")
	os.system("rm *loci*")

	print("{} records had no variant sites or no parsimony informative sites and were removed".format(len(no_variant_no_parsimony_lst)))
	print("{} records remaining in passed_rnas list".format(len(passed_rnas)))
	print("{} records remaining in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))
	print("{} records remaining in cds_parent_rna_dict".format(len(cds_parent_rna_dict)))

def get_cds_lengths():
	get_file_lst = 'find ./vcf2fasta_CDS/ -type f -name "*.fas" > fasta_file_list'
	os.system(get_file_lst)

	cds_length_dict = dict()
	not_multiple_of_three_counter = 0
	not_multiple_of_three_lst = []
	passed_cds_length_dict = dict()

	files = open("fasta_file_list", "r").read().splitlines()

	os.system("rm fasta_file_list")

	for file in files:
		cds = file.split("/")[-1].split(".fas")[0]
		length = len(open(file,"r").read().splitlines()[1])
		cds_length_dict[cds] = length

	for key,value in cds_length_dict.items():
		if value % 3 != 0:
			not_multiple_of_three_counter += 1
			not_multiple_of_three_lst.append(key)
		else:
			passed_cds_length_dict[key] = value

	os.system("mkdir not_multiple_of_three")
	for cds in not_multiple_of_three_lst:
		rna = cds_parent_rna_dict[cds]
		passed_rnas.remove(rna)
		filtered_mrna_gene_dict.pop(rna)
		cds_parent_rna_dict.pop(cds)
		move = "mv vcf2fasta_CDS/{}.fas not_multiple_of_three/".format(cds)
		os.system(move)

	print("Number of concatenated CDS records: {}".format(len(cds_length_dict)))
	print("Number of concatenated CDS records that are not a multiple of 3: {}".format(not_multiple_of_three_counter))
	
	with open("passed_cds_length_stats.txt", "a") as f:
		f.write("Mean cds length: {}".format(str(statistics.mean(passed_cds_length_dict.values()))) + "\n")
		f.write("Median cds length: {}".format(str(statistics.median(passed_cds_length_dict.values()))) + "\n")
		f.write("Minimum cds length: {}".format(str(min(passed_cds_length_dict.values()))) + "\n")
		f.write("Number cds shorter than 2000 base pairs: {}".format(str(len([value for value in passed_cds_length_dict.values() if value <= 2000]))) + "\n")

	with open("passed_cds_length_dist.txt", "a") as f2:
		for value in passed_cds_length_dict.values():
			f2.write(str(value) + "\n")

	with open("passed_cds_lengths.txt", "a") as f3:
		for key,value in passed_cds_length_dict.items():
			f3.write(key + "\t" + str(value) + "\n")

	os.system("cat passed_cds_lengths.txt | sort -k2,2n > passed_CDS_lengths.txt")
	os.system("rm passed_cds_lengths.txt")

# Write csv file of RNAs passing all filters and their coverage metrics
def write_passed_rna_dict_csv():
	csv_file = open("passed_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	records_written = 0 
	
	for i in range(4):
		for sample in subset_sample_list:
			header.append(sample)
	
	writer.writerow(header)

	for key,value in passed_rna_dict.items():
		if key in filtered_mrna_gene_dict.keys():
			row = [key] + value[0] + value[1] + value[2] + value[3]
			records_written += 1
			writer.writerow(row)

	csv_file.close()

	print("{} records written to passed_rna.csv".format(records_written))

	print("{} mRNAs in passed_rnas list".format(len(passed_rnas)))
	print("{} mRNAs in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))

	with open("passed_rnas.txt","a") as f:
		for record in passed_rnas:
			f.write(record + "\n")

def run_iqtree():
	run_iqtree = "iqtree2 -S vcf2fasta_CDS/ -m MFP -T 8"
	os.system(run_iqtree)

def clean_up_iqtree_files():
	delete_gz = 'find ./vcf2fasta_CDS/ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./vcf2fasta_CDS/ -type f -name "*.bionj" -delete'
	delete_contree = 'find ./vcf2fasta_CDS/ -type f -name "*.contree" -delete'
	delete_iqtree = 'find ./vcf2fasta_CDS/ -type f -name "*.iqtree" -delete'
	delete_log = 'find ./vcf2fasta_CDS/ -type f -name "*.log" -delete'
	delete_mldist = 'find ./vcf2fasta_CDS/ -type f -name "*.mldist" -delete'
	#os.system(delete_gz)
	#os.system(delete_bionj)
	#os.system(delete_contree)
	#os.system(delete_iqtree)
	#os.system(delete_log)
	#os.system(delete_mldist)

	#cat_treefiles = 'find ./vcf2fasta_CDS/ -type f -name "*.treefile" -exec cat {} \\; > loci.treefile'
	#cat_boottrees = 'find ./vcf2fasta_CDS/ -type f -name "*.boottrees" -exec cat {} \\; > loci.boottrees'
	
	#os.system(cat_treefiles)
	#os.system(cat_boottrees)

	#delete_treefile = 'find ./vcf2fasta_CDS/ -type f -name "*.treefile" -delete'
	#delete_boottrees = 'find ./vcf2fasta_CDS/ -type f -name "*.boottrees" -delete'
	
	#os.system(delete_treefile)
	#os.system(delete_boottrees)

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
# nw_topology creates cladogram. Option -I gets rid of branc lengths. 
# nw_order orders the tree so that trees with identical topologies will have identical newick strings. Ooption -c d reorders the tree in such a way as to remove the ladder. 
def clean_gene_trees(input_file, output_file, outgroup):
	clean = "{}nw_reroot {} {} | {}nw_topology -I - | {}nw_order -c d - | {}nw_prune - franciscanus > {}".format(nw_utils, input_file, outgroup, nw_utils, nw_utils, nw_utils, output_file)
	os.system(clean)

def main():
	#subset_coverage_dict()

	#bed_file_list = get_zipped_bed_file_list()
	
	#initialize_rna_dict()

	#for regions_file, thresholds_file in bed_file_list:
		#for sample in subset_sample_list:
			#if sample in regions_file and sample in thresholds_file:
				#fill_rna_dict(regions_file, thresholds_file)

	#write_all_rna_dict_csv()
	#filter_rna_dict()
	#get_mrna_gff()
	#create_scaffold_dict()
	#check_proximity()
	#get_passed_rnas()
	#write_new_bed_file()

	#gene_ids = get_gene_ids()

	#os.system("mkdir single_gene_gff_records/")
	#Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	
	# Concatenate all single gene gff records into "sco_gff.gff" file
	#os.system('find ./single_gene_gff_records/ -type f -name "*.record" -exec cat {} \\; > sco_gff.gff')
	
	# Delete the single gene records
	#os.system('find ./single_gene_gff_records/ -type f -name "*.record" -delete')
	#os.system('rmdir single_gene_gff_records/')

	#run_vcf2fasta()

	#remove_redundant_isoforms()

	#replace_missing_genotype_char()
	#identify_no_variant_no_parsimony()
	#remove_no_variant_no_parsimony()
	#get_cds_lengths()
	#write_passed_rna_dict_csv()

	run_iqtree()

	#clean_up_iqtree_files()
	edit_tree_files("loci.treefile","single_locus_trees.nwk")
	#edit_tree_files("loci.boottrees", "single_locus_trees_boot.nwk")
	clean_gene_trees("single_locus_trees.nwk", "clean_trees.nwk", "franciscanus")
	#clean_gene_trees("single_locus_trees_boot.nwk", "clean_trees_boot.nwk", "franciscanus")

if __name__ == "__main__":
	main()
