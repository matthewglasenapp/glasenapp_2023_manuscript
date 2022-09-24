#!/usr/bin/env python

import re

# Create new file that does not have lines beginning with "Warning"
file1 = open('loci.log','r')
file2 = open('test_remove_warning.txt','w')
for line in file1.readlines():
	x = re.findall("^WARNING:", line)
	if not x:
		file2.write(line)

file1.close()
file2.close()

lines = open("test_remove_warning.txt","r").read().splitlines()
gene_list = {}

failed_genes = {}
failed_gene_names = []
passed_genes =[]

for item in lines:
	if item[0:22] == "Reading alignment file":
		index = lines.index(item)
		gene = item.split("/")[1].split(".fas")[0]
		print(gene)
		gene_file = item.split("file ")[1].split(" ...")[0]
		print(gene_file)
		number_sites = lines[index+4].split(" ")[5]
		print(number_sites)
		parsimony_informative_sites = lines[index + 5].split(" ")[0]
		print(parsimony_informative_sites)
		gap_ambiguity_1 = lines[index+9].split("%")[0].split(" ")[-1]
		gap_ambiguity_2 = lines[index+10].split("%")[0].split(" ")[-1]
		gap_ambiguity_3 = lines[index+11].split("%")[0].split(" ")[-1]
		gap_ambiguity_4 = lines[index+12].split("%")[0].split(" ")[-1]
		gap_ambiguity_5 = lines[index+13].split("%")[0].split(" ")[-1]
		gap_ambiguity_6 = lines[index+14].split("%")[0].split(" ")[-1]
		gap_ambiguity_7 = lines[index+15].split("%")[0].split(" ")[-1]
		gap_ambiguity_8 = lines[index+16].split("%")[0].split(" ")[-1]
		gap_ambiguity_9 = lines[index+17].split("%")[0].split(" ")[-1]
		gap_ambiguity_10 = lines[index+18].split("%")[0].split(" ")[-1]

		gene_list[gene] = [number_sites,parsimony_informative_sites,gap_ambiguity_1,gap_ambiguity_2,gap_ambiguity_3,gap_ambiguity_4]

		if int(float(number_sites)) < 300 or int(float(parsimony_informative_sites)) == 0 or int(float(gap_ambiguity_1)) > 50 or int(float(gap_ambiguity_2)) > 50 or int(float(gap_ambiguity_3)) > 50 or int(float(gap_ambiguity_4)) > 50:
			failed_genes[gene] = [number_sites,parsimony_informative_sites,gap_ambiguity_1,gap_ambiguity_2,gap_ambiguity_3,gap_ambiguity_4]
			failed_gene_names.append(gene)
		else:
			passed_genes.append(gene_file)

print(len(gene_list))
print(len(failed_genes))
print(len(failed_gene_names))
print(len(passed_genes))

failed_genes_file = open("failed_genes.txt","a")
for gene in failed_gene_names:
	failed_genes_file.write(gene + "\n")
failed_genes_file.close()

passed_genes_file = open("passed_genes.txt","a")
for gene in passed_genes:
	passed_genes_file.write(gene + "\n")
passed_genes_file.close()



