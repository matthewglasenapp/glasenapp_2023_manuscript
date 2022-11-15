"""
This script calculate coverage depth metrics for mRNA transcripts using mosdepth. 
Mosdepth is used to estimate depth for all protein coding exons.
The exons of each transcript are aggregated to calcuate depths metrics for each mRNA transcript. 
The mosdepth command that was run to generate the input files was:
mosdepth --by regions_file --no-per-base --thresholds 1,5,10,20,30,100 -t 4 --fast-mode prefix bam_file
"""

import gzip 
import os

# Specify directory containing mosdepth outputfiles (both regions.bed.gz and thresholds.bed.gz files)
bed_file_dir = "/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/"

# Initialize dictionary for mRNA names and their average coverage.
rna_dict_depth = dict()

# Initialize dictionary for mRNA names and the number of bases covered by 1, 10, and 20 reads. 
rna_dict_threshold = dict()

# Get zipped list of mosdepth regions and thresholds files for each species
def get_coverage_files_paths():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* | grep -v 'csi' > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	return list(file_list)

# Calculate coverage depth metrics for mRNA transcripts 
def get_mRNA_cov(regions_file, thresholds_file):
	
	line_index = 0
	
	with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file,"rt") as f2:
		
		records = f1.read().splitlines()
		thresholds = f2.read().splitlines()[1:]

		zipped_list = list(zip(records,thresholds))
	
		for record, threshold in zipped_list:

			# Skip over an exon record if it has already been accounted for
			if line_index > records.index(record) and line_index != 0:
				continue
		
			else:
			
				current_mrna = record.split("\t")[3].split("exon-")[1].split("-")[0]
			
				exon_dict = dict()

				current_exon = record.split("\t")[3]
				current_exon_length = float(record.split("\t")[2]) - float(record.split("\t")[1])
				current_exon_coverage = float(record.split("\t")[4])
				
				num_1x = int(threshold.split("\t")[4])
				num_10x = int(threshold.split("\t")[6])
				num_20x = int(threshold.split("\t")[7])

				exon_dict[current_exon] = [current_exon_length, current_exon_coverage, num_1x, num_10x, num_20x]
			
				try:
					next_record = records[line_index + 1]
					next_threshold = thresholds[line_index + 1]
					next_mrna = next_record.split("\t")[3].split("exon-")[1].split("-")[0]
			
				except IndexError:
					rna_dict_depth[current_mrna] = current_exon_coverage
					break
	
				if current_mrna == next_mrna:
				
					while current_mrna == next_mrna:
						# Add additional exons to exons dict 
						# {"Exon_name": [exon_length, exon_coverage, num_1x, num_10x, num_20x]
						exon_dict[next_record.split("\t")[3]] = [(float(next_record.split("\t")[2]) - float(next_record.split("\t")[1])), float(next_record.split("\t")[4]), int(next_threshold.split("\t")[4]), int(next_threshold.split("\t")[6]), int(next_threshold.split("\t")[7])]
				
						line_index += 1
					
						try:
							next_mrna = records[line_index + 1].split("\t")[3].split("exon-")[1].split("-")[0]
							next_threshold = thresholds[line_index + 1]
							next_record = records[line_index + 1]
					
						except IndexError:
							break

					mrna_length = sum([item[0] for item in exon_dict.values()])

					total_1x = sum([item[2] for item in exon_dict.values()])
					total_10x = sum([item[3] for item in exon_dict.values()])
					total_20x = sum([item[4] for item in exon_dict.values()])

					mean = 0 

					for value in exon_dict.values():
						# Calculate weighted average. Sum of (coverage * (exon length/mRNA length))
						mean += (value[1] * (value[0]/mrna_length))

					# {"mRNA_name": mean coverage}
					rna_dict_depth[current_mrna] = mean

					rna_dict_threshold[current_mrna] = [mrna_length, total_1x, total_10x, total_20x]

					line_index +=1 
			
				else:
					rna_dict_depth[current_mrna] = current_exon_coverage

					rna_dict_threshold[current_mrna] = [current_exon_length, num_1x, num_10x, num_20x]
					
					line_index += 1

# Write results to output files 
def write_results(out1, out2):
	with open(out1,"a") as f:
		for key,value in rna_dict_depth.items():
			f.write(str(key) + "\t" + str(value) + "\n")

	with open(out2, "a") as f2:
		for key,value in rna_dict_threshold.items():
			f2.write(str(key) + "\t" + str(int(value[0])) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\n")


def main():
	#coverage_file_list = get_coverage_files_paths()

	# Output from get_coverage_files_paths()
	# get_coverage_files_paths() throws an error when the script is submitted as a slurm array job
	coverage_file_list = [('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/depressus_SRR5767284.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/depressus_SRR5767284.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/droebachiensis_SRR5767286.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/droebachiensis_SRR5767286.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/fragilis_SRR5767279.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/fragilis_SRR5767279.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/franciscanus_SRR5767282.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/franciscanus_SRR5767282.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/intermedius_SRR5767280.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/intermedius_SRR5767280.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/lividus_ERS2351987.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/lividus_ERS2351987.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/nudus_SRR5767281.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/nudus_SRR5767281.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pallidus_SRR5767285.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pallidus_SRR5767285.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pulcherrimus_DRR107784.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pulcherrimus_DRR107784.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pulcherrimus_SRR5767283.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/pulcherrimus_SRR5767283.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/purpuratus_SRR6281818.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/purpuratus_SRR6281818.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/purpuratus_SRR7211988.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/purpuratus_SRR7211988.thresholds.bed.gz'), ('/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/variegatus_SRR7207203.regions.bed.gz', '/hb/home/mglasena/dissertation/data/mosdepth/mosdepth_exons/variegatus_SRR7207203.thresholds.bed.gz')]	
	
	array_id = os.environ["array_id"]
	print("Array ID: {}".format(array_id))

	regions_file = coverage_file_list[int(array_id)][0]
	thresholds_file = coverage_file_list[int(array_id)][1]
	print(regions_file)
	print(thresholds_file)

	print("Caclulating mean coverage for mRNA molecules in {}".format(regions_file))

	get_mRNA_cov(regions_file, thresholds_file)
	
	out1 = regions_file.split("/")[-1].split(".gz")[0]
	out2 = thresholds_file.split("/")[-1].split(".gz")[0]
	write_results(out1, out2)

if __name__ == "__main__":
	main()
