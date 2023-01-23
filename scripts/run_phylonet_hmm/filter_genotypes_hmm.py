import os 

reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Raw genotype calls file
#genotype_calls = "/hb/groups/pogson_group/dissertation/data/raw_vcf_files/genotype_calls.g.vcf.gz"

genotype_calls_split_multiallelics = "/hb/groups/pogson_group/dissertation/data/raw_vcf_files/genotype_calls_split_multiallelics.g.vcf.gz"

output_directory = "/hb/scratch/mglasena/data/genotypes/phylonet_hmm/dro-pal/"

samples_to_include = {
"fragilis_SRR5767279" : "QB3KMK013",
#"nudus_SRR5767281" : "QB3KMK011",
#"franciscanus_SRR5767282" : "QB3KMK010",
#"depressus_SRR5767284" : "QB3KMK015",
"pallidus_SRR5767285" : "QB3KMK002",
"droebachiensis_SRR5767286" : "QB3KMK014",
#"purpuratus_SRR6281818" : "S.purpuratus#1",
#"lividus_ERS2351987" : "4",
"pulcherrimus_SRR5767283" : "QB3KMK016",
#"intermedius_SRR5767280" : "QB3KMK012",
#"purpuratus_SRR7211988" : "SPUR.00",
#"pulcherrimus_DRR107784" : "SAMD00098133"
}

def separate_SNP_INDEL():
	sample_string = ""
	for sample in samples_to_include.values():
		sample_string += "-sn " + sample + " "
	sample_string = sample_string.strip()
	
	input_file = genotype_calls_split_multiallelics
	output_snp = output_directory + "genotype_calls_snv.g.vcf.gz"
	output_indel = output_directory + "genotype_calls_indel.g.vcf.gz"
	get_snp = "gatk SelectVariants -V {} {} --select-type-to-include SNP --output {}".format(input_file, sample_string, output_snp)
	get_indel = "gatk SelectVariants -V {} {} --select-type-to-include INDEL --output {}".format(input_file, sample_string, output_indel)
	os.system(get_snp)
	os.system(get_indel)

def filter_variants():
	input_snp = output_directory + "genotype_calls_snv.g.vcf.gz"
	input_indel = output_directory + "genotype_calls_indel.g.vcf.gz"
	output_snp = output_directory + "filtered_snv.g.vcf.gz"
	output_indel = output_directory + "filtered_indel.g.vcf.gz"

	#Alternate option: bcftools filter -e 'QUAL<30 || FS>60 || SOR>3 || MQ<40 || MQRankSum<-12.5 || QD<2 || ReadPosRankSum<-8' -O z -o output input
	filter_SNPs = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'.format(output_snp, input_snp)

	filter_indels = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'.format(output_indel, input_indel)

	os.system(filter_SNPs)
	os.system("rm " + input_snp)
	os.system(filter_indels)
	os.system("rm " + input_indel)

def merge_vcfs():
	input_snp = output_directory + "filtered_snv.g.vcf.gz"
	input_indel = output_directory + "filtered_indel.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	
	merge_vcfs = "gatk MergeVcfs -I {} -I {} -O {}".format(input_snp, input_indel, output_file)
	os.system(merge_vcfs)
	os.system("rm " + input_snp)
	os.system("rm " + input_indel)

# Set individual genotypes with low quality or read depth to missing: -S . -e 'FMT/DP<3 | FMT/GQ<20'
# Filter SNPs within 3 base pairs of indel: --SnpGap 3
# Remove monomorphic SNPs where no alternative alleles are called for any of the samples: -e 'AC==0'
# Remove insertions from vcf file:
def bcftools_filter():
	input_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	output_file = output_directory + "3bp_filtered_genotype_calls.g.vcf.gz"
	filter = "bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' -Ou {} | bcftools filter --SnpGap 3 -e 'AC==0' -Oz -o {}".format(input_file, output_file)
	os.system(filter)
	os.system("rm " + input_file)

def select_passed_variants():
	filtered_vcf = "3bp_filtered_genotype_calls.g.vcf.gz"
	output_dir = output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz"
	select_variants = "gatk SelectVariants -R {} -V {} -O {} -select-type SNP --exclude-filtered true".format(reference_genome, filtered_vcf, output_dir)
	os.system(select_variants)
	# Alternative option: cat $filtered_vcf | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' > my_PASS.vcf
	#os.system("rm " + filtered_vcf)

def index_vcf(input_file):
	index = "gatk IndexFeatureFile -I {}".format(input_file)
	os.system(index)

def vcf_stats(input_file):
	get_samples_file = "bcftools query -l {} > samples_file.txt".format(input_file)
	os.system(get_samples_file)
	get_stats = "bcftools stats --samples-file samples_file.txt {}".format(input_file)
	os.system(get_stats)
	os.system("rm samples_file.txt")

def main():
	#separate_SNP_INDEL()
	#filter_variants()
	#merge_vcfs()
	bcftools_filter()
	index_vcf(output_directory + "3bp_filtered_genotype_calls.g.vcf.gz")
	select_passed_variants()
	vcf_stats(output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz")

if __name__ == "__main__":
	main()