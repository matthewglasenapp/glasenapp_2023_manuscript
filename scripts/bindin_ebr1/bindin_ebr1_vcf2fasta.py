import os
import subprocess

vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"
vcf_file = "/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls.g.vcf.gz"
gff_file = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"
feature = "gene"

def subset_gff():
	grep_bindin = 'grep "Dbxref=GeneID:373276" {} > bindin.gff'.format(gff_file)
	grep_ebr1 = 'grep "Dbxref=GeneID:577775" {} > ebr1.gff'.format(gff_file)
	os.system(grep_bindin)
	os.system(grep_ebr1)

	concatenate_ebr1_bindin = "cat bindin.gff ebr1.gff > bindin_ebr1.gff"
	os.system(concatenate_ebr1_bindin)

	os.system("rm bindin.gff")
	os.system("rm ebr1.gff")

def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff bindin_ebr1.gff --feat {}".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)
	
def run_iqtree():
	os.system("find vcf2fasta_gene/ -type f > files")
	with open("files","r") as f:
		fasta_alignment_file_list = f.read().splitlines()
	os.system("rm files")

	for file in fasta_alignment_file_list:
		run_iqtree = "iqtree -s {} -m MFP -B 10000 -o QB3KMK010".format(file)
		os.system(run_iqtree)

def main():
	subset_gff()
	run_vcf2fasta()
	replace_missing_genotype_char()
	run_iqtree()

if __name__ == "__main__":
	main()