import gzip 
import os

bed_file_dir = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/"

# Initialize dictionary for mRNA names and their average coverage
rna_dict = dict()

def get_coverage_files_paths():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

def get_mRNA_cov(input_file):
	
	line_index = 0
	
	records = gzip.open(input_file, "rt").read().splitlines()
	
	for record in records:
		
		if line_index > records.index(record) and line_index != 0:
			continue
		
		else:
			
			current_mrna = record.split("\t")[3].split("exon-")[1].split(".")[0]
			
			exon_dict = dict()

			current_exon = record.split("\t")[3]
			current_exon_length = float(record.split("\t")[2]) - float(record.split("\t")[1])
			current_exon_coverage = float(record.split("\t")[4])

			exon_dict[current_exon] = [current_exon_length, current_exon_coverage]
			
			try:
				next_record = records[line_index + 1]
				next_mrna = next_record.split("\t")[3].split("exon-")[1].split(".")[0]
			
			except IndexError:
				rna_dict[current_mrna] = current_exon_coverage
				break
	
			if current_mrna == next_mrna:
				
				while current_mrna == next_mrna:
					# Add additional exons to exons dict 
					# {"Exon_name": [exon_length, exon_coverage]
					exon_dict[next_record.split("\t")[3]] = [float(next_record.split("\t")[2]) - float(next_record.split("\t")[1]), float(next_record.split("\t")[4])]
				
					line_index += 1
					
					try:
						next_mrna = records[line_index + 1].split("\t")[3].split("exon-")[1].split(".")[0]
						next_record = records[line_index + 1]
					
					except IndexError:
						break

				mrna_length = sum([item[0] for item in exon_dict.values()])

				mean = 0 

				for value in exon_dict.values():
					# Calculate weighted average. Sum of (coverage * exon length)/mRNA length)
					mean += (value[1] * (value[0]/mrna_length))

				# {"mRNA_name": mean coverage}
				rna_dict[current_mrna] = mean

				line_index +=1 
			
			else:
				rna_dict[current_mrna] = current_exon_coverage
				line_index += 1

def write_results(output_file):
	with open(output_file,"a") as f:
		for key,value in rna_dict.items():
			f.write(str(key) + "\t" + str(value) + "\n")

def main():
	coverage_file_list = open("regions_files","r").read().splitlines()

	array_id = os.environ["array_id"]
	print("Array ID: {}".format(array_id))
	
	input_file = coverage_file_list[int(array_id)]
	print("Caclulating mean coverage for mRNA molecules in {}".format(input_file))

	get_mRNA_cov(input_file)
	
	output_file = input_file.split(".regions")[0] + "_mrna_cov.tsv"
	write_results(output_file)

if __name__ == "__main__":
	main()
