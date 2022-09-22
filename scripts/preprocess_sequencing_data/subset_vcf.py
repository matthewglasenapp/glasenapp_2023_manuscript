import os 

filtered_vcf = "/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls.g.vcf.gz"

samples_to_include = {
"fragilis_SRR5767279" : "QB3KMK013",
#"nudus_SRR5767281" : "QB3KMK011",
"franciscanus_SRR5767282" : "QB3KMK010",
#"depressus_SRR5767284" : "QB3KMK015",
"pallidus_SRR5767285" : "QB3KMK002",
"droebachiensis_SRR5767286" : "QB3KMK014",
#"purpuratus_SRR6281818" : "S.purpuratus#1",
#"lividus_ERS2351987" : "4",
"pulcherrimus_SRR5767283" : "QB3KMK016",
"intermedius_SRR5767280" : "QB3KMK012",
"purpuratus_SRR7211988" : "SPUR.00",
#"pulcherrimus_DRR107784" : "SAMD00098133"
}

def subset_vcf():
	sample_string = ""
	for sample in samples_to_include.values():
		sample_string += "-sn " + sample + " "
	sample_string = sample_string.strip()
	
	input_file = filtered_vcf
	output_file = "subset_3bp_filtered_genotype_calls.g.vcf.gz"
	subset = "gatk SelectVariants -V {} {} --output {}".format(input_file, sample_string, output_file)
	os.system(subset)

def main():
	subset_vcf()

if __name__ == "__main__":
	main()