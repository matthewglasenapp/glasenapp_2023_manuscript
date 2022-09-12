import os 

root_dir = "/hb/scratch/mglasena/data/"

# Path to S. purpuratus reference genome file
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Directory for single sample vcf files
vcf_dir = root_dir + "vcf_files/"

# Directory for multisample vcf 
multisample_vcf_dir = root_dir + "combined_vcf/"
make_multisample_vcf_dir = "mkdir -p {}".format(multisample_vcf_dir)
os.system(make_multisample_vcf_dir)

def combine_GVCFs():
	get_file_paths = "find {} -type f -name '*.gz' | grep -v 'tbi' > single_sample_vcf_file_paths.txt".format(vcf_dir)
	os.system(get_file_paths)
	file_paths_string = " ".join(["-V " + path for path in open("single_sample_vcf_file_paths.txt","r").read().splitlines()])
	os.system("rm single_sample_vcf_file_paths.txt")
	output_file = multisample_vcf_dir + "raw_combined_vcf.g.vcf.gz"
	combine_gvcfs = "gatk CombineGVCFs -R {} -O {} {}".format(reference_genome, output_file, file_paths_string)
	os.system(combine_gvcfs)

def genotype_GVCFs():
	input_file = multisample_vcf_dir + "raw_combined_vcf.g.vcf.gz"
	output_file = multisample_vcf_dir + "genotype_calls.g.vcf.gz"
	call_genotypes = "gatk GenotypeGVCFs -R {} -V {} -O {}".format(reference_genome, input_file, output_file)
	os.system(call_genotypes)

def main():
	combine_GVCFs()
	genotype_GVCFs()

if __name__ == "__main__":
	main()