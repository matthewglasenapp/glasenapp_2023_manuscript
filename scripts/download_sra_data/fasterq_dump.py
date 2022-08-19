"""
Download multiple accession from the NCBI SRA. Consider modifying to parallelize.
"""

# SRR5767283 (pulcherrimus) and SRR5767280 (intermedius) are skipped because they were not uploaded to ncbi correctly. I have the raw sequencing reads for these two samples. 

import os

prefetch = "/hb/groups/pogson_group/dissertation/software/sratoolkit.2.11.2-centos_linux64/bin/prefetch"
fasterq_dump = "/hb/groups/pogson_group/dissertation/software/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump"

threads = 8

prefetch_directory = "/hb/scratch/mglasena/prefetch/"
make_prefetch_dir = "mkdir -p {}".format(prefetch_directory)
os.system(make_prefetch_dir)

output_directory = "/hb/home/mglasena/short_read_data/"
make_output_dir = "mkdir -p {}".format(output_directory)
os.system(make_output_dir)

temporary_directory = "/hb/scratch/mglasena/SRA/"

# skipped_accessions = {"SRR5767283": "pulcherrimus,QB3KMK016,SAMN07269099", "SRR5767280": "intermedius,QB3KMK012,SAMN07269102"}

# species, sample_name, ncbi biosample

#urchin_sra_accessions = {
	#"SRR5767279": "fragilis,QB3KMK013,SAMN07269103",
	#"SRR5767281": "nudus,QB3KMK011,SAMN07269101",
	#"SRR5767282": "franciscanus,QB3KMK010,SAMN07269100",
	#"DRR107784": "pulcherrimus,SAMD00098133,SAMD00098133",
	#"SRR5767284": "depressus,QB3KMK015,SAMN07269098",
	#"SRR5767285": "pallidus,QB3KMK002,SAMN07269097",
	#"SRR5767286": "droebachiensis,QB3KMK014,SAMN07269096",
	#"SRR6281818": "purpuratus,S.purpuratus#1,SAMN08013506",
	#"SRR7211988": "purpuratus,SPUR.00,SAMN00829422",
	#"ERR5621404": "lividus,4,ERS2351987"
	#"SRR2583947": "franciscanus,Sf1,SAMN04156173"}

urchin_sra_accessions = {
	#"ERR5621404": "lividus", "4", "ERS2351987",
	"ERR5621405": "lividus,4,ERS2351987",
	"ERR5621406": "lividus,4,ERS2351987",
	"ERR5621407": "lividus,4,ERS2351987",
	"ERR5621408": "lividus,4,ERS2351987",
	"ERR5621409": "lividus,4,ERS2351987",
	"ERR5621410": "lividus,4,ERS2351987",
	"ERR5621411": "lividus,4,ERS2351987",
	"ERR5621412": "lividus,4,ERS2351987",
	"ERR5621413": "lividus,4,ERS2351987",
	"ERR5621414": "lividus,4,ERS2351987",
	"ERR5621415": "lividus,4,ERS2351987",
	"ERR5621416": "lividus,4,ERS2351987",
	"ERR5621417": "lividus,4,ERS2351987",
	"ERR5621418": "lividus,4,ERS2351987",
}

def prefetch(accession_number):
	output_dir = prefetch_directory + accession_number
	"{} {} -O {} --max-size u".format(prefetch, accession_number, output_dir)

def download_accession(accession_number):
	prefetch_dir = prefetch_directory + accession_number
	os.chdir(prefetch_dir)
	temp_dir = temporary_directory + accession_number 
	download = "{} {} --outdir {} --split-3 -t {} -e {}".format(fasterq_dump, accession_number, output_directory, temp_dir, threads)

	os.system(download)

	os.chdir(output_directory)

	# Remove prefetch dir
	remove_prefetch = "rm -r {}".format(prefetch_dir)
	os.sytem(remove_prefetch)

	# Remove temp_dir
	remove_temp = "rm -r {}".format(temp_dir)
	os.system(remove_temp)

	output_file_1 = accession_number + "_1.fastq"
	output_file_2 = accession_number + "_2.fastq"
	gzip_1 = "gzip {}".format(output_file_1)
	gzip_2 = "gzip {}".format(output_file_1)
	os.system(gzip_1)
	os.system(gzip_2)	

def main():
	for key in urchin_sra_accessions:
		print("Prefetching {}".format(key))
		prefetch(key)
		print("{} prefetch complete".format(key))
		
		print("Downloading {}".format(key))
		download_accession(key)
		print("{} complete".format(key))

if __name__ == "__main__":
	main()
