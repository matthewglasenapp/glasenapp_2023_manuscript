#!/usr/bin/env python

import re

iqtree_log = "loci.log"

number_samples_analyzed = 10

gene_dict = {}
passed_genes = {}
passed_gene_names = []
passed_gene_files = []
failed_genes = {}
failed_gene_names = []
failed_gene_files = []

# Create new file that does not have lines beginning with "Warning"
def remove_warning_lines():
	with open(iqtree_log,'r') as file1:
		with open('warning_removed.txt','w') as file2:
			for line in file1.readlines():
				x = re.findall("^WARNING:", line)
				if not x:
					file2.write(line)

def parse_log_file():
	lines = open("warning_removed.txt","r").read().splitlines()
	for line in lines:
		if line[0:22] == "Reading alignment file":
			index = lines.index(line)
			gene_name = line.split("/")[1].split(".fas")[0]
			gene_file = line.split("file ")[1].split(" ...")[0]
			number_sites = int(lines[index+4].split(" ")[5])
			parsimony_informative_sites = int(lines[index + 5].split(" ")[0])

			gap_ambiguity_list = []
			for n in range(index+8, (index + 8 + number_samples_analyzed)):
				gap_ambiguity_percent = float(lines[n].split("%")[0].split(" ")[-1])
				gap_ambiguity_list.append(gap_ambiguity_percent)

			# Add data for all genes to gene_dict
			gene_dict[gene_name] = [number_sites,parsimony_informative_sites,gap_ambiguity_list]

			# If no parsimony informative sites or gap_ambiguity > 50, add to failed_genes dict
			#if parsimony_informative_sites == 0 or max(gap_ambiguity_list) > 50:
			if max(gap_ambiguity_list) > 50:
				failed_genes[gene_name] = [number_sites,parsimony_informative_sites,gap_ambiguity_list]
				failed_gene_names.append(gene_name)
				failed_gene_files.append(gene_file)
			else:
				passed_genes[gene_name] = [number_sites,parsimony_informative_sites,gap_ambiguity_list]
				passed_gene_files.append(gene_file)
				passed_gene_names.append(gene_name)

def write_failed_genes_file():
	with open("failed_genes.txt","a") as failed_genes_file:
		for gene in failed_gene_files:
			failed_genes_file.write(gene + "\n")

def write_passed_genes_file():
	with open("passed_genes.txt","a") as passed_genes_file:
		for gene in passed_gene_files:
			passed_genes_file.write(gene + "\n")

def main():
	remove_warning_lines()
	parse_log_file()
	write_failed_genes_file()
	write_passed_genes_file()

if __name__ == "__main__":
	main()