"""
This script takes a fasta file of unaligned sequences, aligns them, and creates a maximum-likelihood tree. 
Speicfy raw unaligned sequence fasta file and number of threads in "input_fasta_file" and "threads"

1) Reformat fasta files downloaded from NCBI. Change from >accession to >species_accession for phylogenetic tree visualization. 
2) Align reformatted sequences using clustal omega 1.2.3
3) Run iqtree on aligned sequences
"""
import os 

# Raw fasta input file 
input_fasta_file = "raw_mitochondrial_sequences.fasta"
threads = 6

# Create new header line for fasta sequence 
def get_new_accession_string(line):
	accession = line.split()[0].split(">")[1]
	species_name = line.split()[2]
	new_accession_name = accession + "_" + species_name
	separator = " "
	new_line = ">" + new_accession_name + " " + separator.join(line.split()[1:])
	return new_line

# Write new fasta file with modified headers
def reformat_fasta(input_file):
	output_file = "reformatted_" + input_fasta_file
	new_file = open(output_file,"a")
	for line in open(input_file,"r").readlines():
		if line[0] == ">":	
			new_file.write(get_new_accession_string(line))
		else:
			new_file.write(line)
	new_file.close()

# Align fasta sequence records
def align_fasta_sequences():
	input_file = "reformatted_" + input_fasta_file
	output_file = "aligned_reformatted_" + input_fasta_file
	align = "clustalo -i {} -o {} --threads={} -v".format(input_file, output_file, threads)
	os.system(align)

# Create phylogenetic tree from aligned sequences
def create_tree():
	input_file = "aligned_reformatted_" + input_fasta_file
	threads = "AUTO"
	bootstrap_replicates = 10000
	run_iqtree = "iqtree -s {} -m MFP -B {} -T {}".format(input_file, bootstrap_replicates, threads)
	os.system(run_iqtree)

def main():
	reformat_fasta(input_fasta_file)
	align_fasta_sequences()
	create_tree()

if __name__ == "__main__":
	main()
