import os 
import subprocess
#import multiprocessing
#from joblib import Parallel, delayed

#num_cores = multiprocessing.cpu_count()

# Get separate multi-sample vcf file for each scaffold. Get coordinates for each site applied to the matrix. 

reference_genome = ""
scaffold_list_file = ""
filtered_vcf = ""
vcf2phylip_path = "/hb/home/mglasena/software/vcf2phylip/"

output_dir = "/hb/groups/pogson_group/vcf2phylip_4way/nexus/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

outgroup_sample_name = "QB3KMK012"
number_species = "4"

def get_scaffold_list():
	f = open(scaffold_list_file, "r")
	scaffold_list = f.read().splitlines()

	return scaffold_list


def subset_vcf_by_scaffold(scaffold):
	subset_vcf_by_scaffold = "gatk SelectVariants -R {} -V {} --select-type-to-include SNP --select-type-to-exclude MNP --select-type-to-exclude INDEL --select-type-to-exclude SYMBOLIC --select-type-to-exclude MIXED -O {}{}/{}.vcf.gz -L {} --select-type-to-exclude NO_VARIATION --exclude-filtered true --exclude-non-variants true --exclude-sample-name QB3KMK011 --exclude-sample-name QB3KMK010 --exclude-sample-name QB3KMK015 --exclude-sample-name QB3KMK016 --exclude-sample-name SPUR.00".format(reference_genome, filtered_vcf, output_dir, scaffold, scaffold, scaffold)
	
	os.system(subset_vcf_by_scaffold)

def convert_vcf_to_nexus(scaffold):
	run_vcf2phylip = "python3 {}vcf2phylip.py -w --input {}{}/{}.vcf.gz --min-samples-locus {} -p -n -r --outgroup {} --output-folder {}{}--output-prefix {}".format(vcf2phylip_path, output_dir, scaffold, scaffold, number_species, outgroup_sample_name, output_dir, scaffold, scaffold)
	
	os.system(run_vcf2phylip)

def reformat_coordinate_files(scaffold):
	os.chdir(output_dir + scaffold)
	
	
	os.environ['coordinate_file'] = "{}.min{}.used_sites.tsv".format(scaffold, number_species)


	create_local_coordinate_file = "cat {}.min{}.used_sites.tsv | grep -v POS | awk '{ print $2 }' > {}_coordinates".format(scaffold, number_species, scaffold)

	create_global_coordinate_file = '''cat $coordinate_file | grep -v POS | awk '{ print $1 ":" $2 }' > coordinates'''
	
	os.system(create_local_coordinate_file)
	subprocess.call(create_global_coordinate_file, shell=True);

def main():
	scaffold_list = get_scaffold_list()
	for scaffold in scaffold_list:
		subset_vcf_by_scaffold(scaffold)
		convert_vcf_to_nexus(scaffold)
		reformat_coordinate_files(scaffold)

if __name__ == "__main__":
        main()
