import os 

# Do I need to normalize after the first time?
#bcftools norm -f /hb/home/mglasena/reference/sp5_0_GCF_genomic.fa -m- -Oz -o norm_joint_genotype.vcf.gz joint_genotype.vcf.gz

threads = 20

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

def bcftools_3bp_filter():
	input_file = multisample_vcf_dir + "genotype_calls.g.vcf.gz"
	output_file = multisample_vcf_dir + "3bp_filter_genotype_calls.g.vcf.gz"
	filter = "bcftools filter -g 3 -O z --output {} {} --threads {}".format(output_file, input_file, threads)
	os.system(bcftools_filter)

def separate_SNP_INDEL():
	input_file = multisample_vcf_dir + "3bp_filter_genotype_calls.g.vcf.gz"
	output_snp = multisample_vcf_dir + "3bp_filter_snv.g.vcf.gz"
	output_indel = multisample_vcf_dir + "3bp_filter_indel.g.vcf.gz"
	get_snp = "gatk SelectVariants -V {} --select-type-to-include SNP --output {}".format(input_file, output_snp)
	get_indel = "gatk SelectVariants -V {} --select-type-to-include INDEL --output {}".format(input_file, output_indel)
	os.system(get_snp)
	os.system(get_indel)
	os.system("rm " + input_file)

def filter_variants():
	input_snp = multisample_vcf_dir + "3bp_filter_snv.g.vcf.gz"
	input_indel = multisample_vcf_dir + "3bp_filter_indel.g.vcf.gz"
	output_snp = multisample_vcf_dir + "filtered_snv.g.vcf.gz"
	output_indel = multisample_vcf_dir + "filtered_indel.g.vcf.gz"

	filter_SNPs = 'gatk VariantFiltration --output {} --variant {} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'.format(output_snp, input_snp)

	filter_indels = 'gatk VariantFiltration --output {} --variant {} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'.format(output_indel, input_indel)

	os.system(filer_SNPs)
	os.system(filter_indels)
	os.system("rm " + input_snp)
	os.system("rm " + input_indel)

def merge_vcfs():
	input_snp = multisample_vcf_dir + "filtered_snv.g.vcf.gz"
	input_indel = multisample_vcf_dir + "filtered_indel.g.vcf.gz"
	output_file = multisample_vcf_dir + "filtered_genotype_calls.g.vcf.gz"
	
	merge_vcfs = "gatk MergeVcfs -I {} -I {} -O {}".format(input_snp, input_indel, output_file)
	os.system(merge_vcfs)
	os.system("rm " + input_snp)
	os.system("rm " + input_indel)

def index_vcf(input_file):
	index = "gatk IndexFeatureFile -I {}".format(input_file)
	os.system(index)

def vcf_stats(input_file):
	get_samples_file = "bcftools query -l {} > samples_file.txt".format(input_file)
	os.system(get_samples_file)
	get_stats = "bcftools stats --samples-file samples_file.txt {}".format(input_file)
	os.system("rm samples_file.txt")

def main():
	#combine_GVCFs()
	genotype_GVCFs()
	index_vcf(multisample_vcf_dir + "genotype_calls.g.vcf.gz")
	#vcf_stats(multisample_vcf_dir + "genotype_calls.g.vcf.gz")
	#bcftools_3bp_filter()
	#separate_SNP_INDEL()
	#filter_variants()
	#merge_vcfs()
	#index_vcf(multisample_vcf_dir + "filtered_genotype_calls.g.vcf.gz")
	#vcf_stats(multisample_vcf_dir + "filtered_genotype_calls.g.vcf.gz")

if __name__ == "__main__":
	main()