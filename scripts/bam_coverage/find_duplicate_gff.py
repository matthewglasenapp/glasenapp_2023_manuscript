import os
import subprocess

def subset_gff():
	subset_gff = '''cat sco_gff | awk '$3 == "gene"' > test'''

	subprocess.call(subset_gff, shell=True);

def get_gene_ids():
	with open("test","r") as f:
		inputs = f.read().splitlines()

	new_inputs = [line.split("\t")[8].split(";")[1] for line in inputs]

	with open("test_2","a") as f2:
		for line in new_inputs:
			f2.write(line + "\n")

	os.system("rm test")

def find_duplicate_ids():
	os.system("sort test_2 | uniq -d > duplicate_records")

	with open("duplicate_records", "r") as f3:
		duplicate_records = f3.read().splitlines()

	print("Duplicate Records:" + "\n")
	
	for record in duplicate_records:
		print(record)

	os.system("rm test_2")
	os.system("rm duplicate_records")

def main():
	subset_gff()
	get_gene_ids()
	find_duplicate_ids()

if __name__ == "__main__":
	main()