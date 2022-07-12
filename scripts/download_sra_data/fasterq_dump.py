"""
Download multiple accession from the NCBI SRA. Consider modifying to parallelize.
"""

# species, sample_name, ncbi biosample, library_name, read_group_string (@<instrument>:<run number>:<flowcell ID>:<lane>)
# For SRA accessions without read group information, the read group string is simply the SRA accession
# For SRA accession without library information, the library is simply the SRA accession
# SRR5767283 (pulcherrimus) and SRR5767280 (intermedius) are skipped because they were not uploaded to ncbi correctly. I have the raw sequencing reads for these two samples. 

import os

skipped_accessions = {"SRR5767283": "pulcherrimus,QB3KMK016,SAMN07269099", "SRR5767280": "intermedius,QB3KMK012,SAMN07269102"}

urchin_sra_accessions = {
	"SRR5767279": "fragilis,QB3KMK013,SAMN07269103",
	"SRR5767281": "nudus,QB3KMK011,SAMN07269101",
	"SRR5767282": "franciscanus,QB3KMK010,SAMN07269100",
	"DRR107786": "pulcherrimus,SAMD00098133,SAMD00098133",
	"SRR5767284": "depressus,QB3KMK015,SAMN07269098",
	"SRR5767285": "pallidus,QB3KMK002,SAMN07269097",
	"SRR5767286": "droebachiensis,QB3KMK014,SAMN07269096",
	"SRR6281818": "purpuratus,S.purpuratus#1,SAMN08013506",
	"SRR7211988": "purpuratus,SPUR.00,SAMN00829422",
	"ERR5671699": "lividus,4,ERS2351987"
	"SRR2583947": "franciscanus,Sf1,SAMN04156173"
}

def download_accession(accession_number):
	os.system("/hb/groups/pogson_group/dissertation/software/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --outdir /hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/ --split-3 -t /hb/scratch/mglasena/ -e 48 " + accession_number)

def main():
	for key in urchin_sra_accessions:
		print("downloading {}".format(key))
		download_accession(key)
		print("{} complete".format(key))

if __name__ == "__main__":
	main()
