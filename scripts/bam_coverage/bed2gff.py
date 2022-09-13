import os 
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

# The script requires genes_pass_filter.bed, the output of running find_sco.py

def get_gene_ids():
	split_columns = "awk '{ print $10 }' genes_pass_filter.bed > gene_list"
	os.system(split_columns)

	with open("gene_list","r") as f:
		with open("gene_ids","a") as f2:
			gene_list = f.read().splitlines()
			for gene in gene_list:
				identifier = gene.split(";")[1]
				f2.write(identifier + "\n")

# time: 23 minutes
#def make_sco_gff():
	#with open("purpuratus.gff", "r") as f:
		#with open("sco_gff", "a") as f2:
			#for line in f:
				#for gene_id in gene_ids:
					#if gene_id in line:
						#f2.write(line)

### Alternative option
# time: 3 minutes 
def make_sco_gff(gene):
	command = "grep {} purpuratus.gff > {}.txt".format(gene, gene)
	os.system(command)

def main():
	get_gene_ids()
	gene_ids = open("gene_ids", "r").read().splitlines()
	#make_sco_gff()
	
	Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	os.system("cat *.txt > sco_gff")
	os.system("rm *.txt")

	# Run command to verify the right number of records have been written
	# cat sco_gff | awk '$3 == "gene"' | wc -l

if __name__ == "__main__":
	main()