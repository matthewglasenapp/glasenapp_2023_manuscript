"""
Download multiple accession from the NCBI SRA. Consider modifying to parallelize.
"""

# SRR5767283 (pulcherrimus) and SRR5767280 (intermedius) are skipped because they were not uploaded to ncbi correctly. I have the raw sequencing reads for these two samples. 

import os

fasterq_dump = "/hb/groups/pogson_group/dissertation/software/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump"

threads = 8

output_directory = "/hb/scratch/mglasena/short_read_data/"
make_output_dir = "mkdir -p {}".format(output_directory)
os.system(make_output_dir)

temporary_directory = "/hb/scratch/mglasena/"

skipped_accessions = {"SRR5767283": "pulcherrimus,QB3KMK016,SAMN07269099", "SRR5767280": "intermedius,QB3KMK012,SAMN07269102"}

# species, sample_name, ncbi biosample

urchin_sra_accessions = {
	#"SRR5767279": "fragilis,QB3KMK013,SAMN07269103",
	#"SRR5767281": "nudus,QB3KMK011,SAMN07269101",
	#"SRR5767282": "franciscanus,QB3KMK010,SAMN07269100",
	"DRR107784": "pulcherrimus,SAMD00098133,SAMD00098133",
	#"SRR5767284": "depressus,QB3KMK015,SAMN07269098",
	#"SRR5767285": "pallidus,QB3KMK002,SAMN07269097",
	#"SRR5767286": "droebachiensis,QB3KMK014,SAMN07269096",
	#"SRR6281818": "purpuratus,S.purpuratus#1,SAMN08013506",
	#"SRR7211988": "purpuratus,SPUR.00,SAMN00829422",
	#"ERR5621404": "lividus,4,ERS2351987"
	#"SRR2583947": "franciscanus,Sf1,SAMN04156173"
}

def download_accession(accession_number):
	temp_dir = temporary_directory + accession_number 
	download = "{} --outdir {} --split-3 -t {} -e {} {}".format(fasterq_dump, output_directory, temp_dir, threads, accession_number)

	os.system(download)

def main():
	for key in urchin_sra_accessions:
		print("downloading {}".format(key))
		download_accession(key)
		print("{} complete".format(key))

if __name__ == "__main__":
	main()
