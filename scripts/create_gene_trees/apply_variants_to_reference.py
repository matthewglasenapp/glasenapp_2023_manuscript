import os

vcf_dir = "/hb/home/mglasena/sample_data/"
input_vcf = "/hb/home/mglasena/sample_data/norm_filtered_calls.vcf.gz"
reference = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

#cat /hb/home/mglasena/reference/sp5_0_GCF_genomic.fa | \
#bcftools consensus \
#--iupac-codes \
#--include 'FILTER="PASS" && strlen(REF)>=strlen(ALT)' \
#--mark-del - \
#-M "N" \
#-H I \
#-s QB3KMK014 \
#-o /hb/home/mglasena/consensus/droebachiensis.fa \
#/hb/home/mglasena/variant/merge/filter/filtered_calls.vcf.gz

dict = {
"SRR5767279" : ["SRR5767279","fragilis","QB3KMK013","SAMN07269103","VJCQB3KMK013","HS3:147:d0gnlacxx:3","d0gnlacxx:3"],
"SRR5767281" : ["SRR5767281","nudus","QB3KMK011","SAMN07269101","VJCQB3KMK011","HS2:148:C0EN2ACXX:4","C0EN2ACXX:4"],
"SRR5767282" : ["SRR5767282","franciscanus","QB3KMK010","SAMN07269100","VJCQB3KMK010","HS2:148:C0EN2ACXX:5","C0EN2ACXX:5"],
#"DRR107786" : ["DRR107786","pulcherrimus","SAMD00098133","SAMD00098133","DRR107786","HWI-ST462R:262:C1J2AACXX:3","C1J2AACXX:3"],
"SRR5767284" : ["SRR5767284","depressus","QB3KMK015","SAMN07269098","VJCQB3KMK015","HS3:171:d0le4acxx:2","d0le4acxx:2"],
"SRR5767285" : ["SRR5767285","pallidus","QB3KMK002","SAMN07269097","VJCQB3KMK002","HS1_0066:8","HS1_0066:8"],
"SRR5767286" : ["SRR5767286","droebachiensis","QB3KMK014","SAMN07269096","VJCQB3KMK014","HS3:171:d0le4acxx:1","d0le4acxx:1"],
#"SRR6281818" : ["SRR6281818","purpuratus","S.purpuratus#1","SAMN08013506","A630_1","SRR6281818","SRR6281818"],
#"ERR5671699" : ["ERR5671699","lividus","4","ERS2351987","4089_HiC_PI_3","ERR5671699","ERR5671699"],
"SRR5767283" : ["SRR5767283","pulcherrimus","QB3KMK016","SAMN07269099","VJCQB3KMK016","HS3:171:d0le4acxx:3","d0le4acxx:3"],
"SRR5767280" : ["SRR5767280","intermedius","QB3KMK012","SAMN07269102","VJCQB3KMK012","HS2:148:C0EN2ACXX:3","C0EN2ACXX:3"],
"SRR7211988" : ["SRR7211988","purpuratus","SPUR.00","SAMN00829422","CIT_GEC_SP_1",["HISEQ:348:H2YWCBCXX:1","HISEQ:348:H2YWCBCXX:2"],["H2YWCBCXX:1","H2YWCBCXX:2"]]
}

def make_consensus(sample):
	output_file = "{}{}.fa".format(vcf_dir, sample)
	consensus = "cat {} | bcftools consensus --iupac-codes --include 'FILTER=\"PASS\" && strlen(REF)>=strlen(ALT)' --mark-del - -M \"N\" -H I -s {} -o {} {}".format(reference, sample, output_file, input_file)
	os.system(consensus)

def main():
	array_id = os.environ["array_id"]
	accession_list = list(dict)
	accession_id = accession_list[int(array_id)]
	sample = dict[accession_id][2]

	make_consensus(input_vcf, sample)

if __name__ == "__main__":
	main()
