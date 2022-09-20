import os
import subprocess

subset_gff = '''cat sco_gff | awk '$3 == "gene"' > test'''

subprocess.call(subset_gff, shell=True);

with open("test","r") as f:
	inputs = f.read().splitlines()

new_inputs = [line.split("\t")[8].split(";")[1] for line in inputs]

with open("test_2","a") as f2:
	for line in new_inputs:
		f2.write(line + "\n")

os.system("sort test_2 | uniq -d")