import gzip 

file = "purpuratus_SRR7211988.regions.bed.gz"
#file = "test.tsv.gz"

def get_mRNA_cov(file):
	rna_dict = dict()
	line_index = 0
	records = gzip.open(file, "rt").read().splitlines()
	for record in records:
		if line_index > records.index(record) and line_index != 0:
			continue
		else:
			current_mrna = record.split("\t")[3].split("exon-")[1].split(".")[0]
			current_exon = record.split("\t")[3]
			current_exon_length = float(record.split("\t")[2]) - float(record.split("\t")[1])
			current_exon_coverage = float(record.split("\t")[4])

			exon_dict = dict()
			exon_dict[current_exon] = [current_exon_length, current_exon_coverage]
			
			try:
				next_record = records[line_index + 1]
				next_mrna = next_record.split("\t")[3].split("exon-")[1].split(".")[0]
			except IndexError:
				rna_dict[current_mrna] = current_exon_coverage
				break
	
			if current_mrna == next_mrna:
				while current_mrna == next_mrna:
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
					mean += value[1] * (value[0]/mrna_length)

				rna_dict[current_mrna] = mean

				line_index +=1 
			
			else:
				rna_dict[current_mrna] = current_exon_coverage
				line_index += 1

	with open("results_mrna.tsv","a") as f:
		for key,value in rna_dict.items():
			f.write(str(key) + "\t" + str(value) + "\n")

def main():
	get_mRNA_cov(file)

if __name__ == "__main__":
	main()
