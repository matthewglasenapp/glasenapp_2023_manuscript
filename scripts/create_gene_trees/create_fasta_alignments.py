import os
import subprocess

vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"
vcf_file = "/hb/scratch/mglasena/data/genotypes/franciscanus/3bp_filtered_genotype_calls.g.vcf.gz"
gff_file = "/hb/scratch/mglasena/test/sco_gff.gff"
feature = "gene"

sample_names = {
'QB3KMK013': 'fragilis',
#'QB3KMK011': 'nudus',
'QB3KMK010': 'franciscanus',
#'QB3KMK015': 'depressus',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
#'S.purpuratus#1': 'purpuratus_SRR6281818',
'QB3KMK016': 'pulcherrimus_SRR5767283',
'QB3KMK012': 'intermedius',
'SPUR.00': 'purpuratus_SRR7211988',
#'SAMD00098133': 'pulcherrimus_DRR107784'
}

def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff {} --feat {}".format(vcf2fasta, reference_genome, vcf_file, gff_file, feature)
	os.system(run_vcf2fasta)

def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)
	
def run_iqtree():
	run_iqtree = "iqtree -S vcf2fasta_gene/ -m MFP --prefix loci"
	os.system(run_iqtree)

def edit_tree_files():
	with open("loci.treefile", "r") as f:
		tree_list = f.read().splitlines()
	
	with open("single_locus_trees.nwk","a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

def main():
	run_vcf2fasta()
	replace_missing_genotype_char()
	run_iqtree()
	edit_tree_files()

if __name__ == "__main__":
	main()