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
bed_file_dir = "/hb/scratch/mglasena/mosdepth/"

mrna_dict = dict()

mrna_dict_depth = dict()

mrna_dict_threshold = dict()

# Get zipped list of mosdepth regions and thresholds files for each species
def get_coverage_files_paths():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* | grep -v 'csi' > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* | grep -v 'csi' > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	os.system("rm regions_files")
	os.system("rm thresholds_files")
	return list(file_list)

def get_mRNA_cov(regions_file, thresholds_file):
	with gzip.open(regions_file, "rt") as f1, gzip.open(thresholds_file,"rt") as f2:
		records = f1.read().splitlines()
		thresholds = f2.read().splitlines()[1:]

		zipped_list = list(zip(records,thresholds))

		for record in zipped_list:
			(mrna, exon, length, coverage, num_1x, num_10x, num_20x) = (record[0].split("\t")[3].split("exon-")[1].split("-")[0], 
				record[0].split("\t")[3], 
				(float(record[0].split("\t")[2]) - float(record[0].split("\t")[1])), 
				float(record[0].split("\t")[4]), 
				int(record[1].split("\t")[4]), 
				int(record[1].split("\t")[5]), 
				int(record[1].split("\t")[6]))

			if mrna in mrna_dict:
				mrna_dict[mrna].append([exon, length, coverage, num_1x, num_10x, num_20x])
			else:
				mrna_dict[mrna] = [[exon, length, coverage, num_1x, num_10x, num_20x]]

	for key,value in mrna_dict.items():
		mrna_length = sum([list[1] for list in value])
		mean_depth = sum([(list[2] * (list[1] / mrna_length)) for list in value])
		total_1x = sum([list[3] for list in value])
		total_10x = sum([list[3] for list in value])
		total_20x = sum([list[5] for list in value])

		mrna_dict_depth[key] = mean_depth

		mrna_dict_threshold[key] = [mrna_length, total_1x, total_10x, total_20x]

def write_results(out1, out2):
	with open(out1,"a") as f:
		for key,value in mrna_dict_depth.items():
			f.write(str(key) + "\t" + str(value) + "\n")

	with open(out2, "a") as f2:
		for key,value in mrna_dict_threshold.items():
			f2.write(str(key) + "\t" + str(int(value[0])) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\n")

def main():
	#print(get_coverage_files_paths())

	coverage_file_list = [('/hb/scratch/mglasena/mosdepth/depressus_SRR5767284.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/depressus_SRR5767284.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/droebachiensis_SRR5767286.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/droebachiensis_SRR5767286.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/fragilis_SRR5767279.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/fragilis_SRR5767279.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/franciscanus_SRR5767282.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/franciscanus_SRR5767282.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/intermedius_SRR5767280.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/intermedius_SRR5767280.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/lividus_ERS2351987.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/lividus_ERS2351987.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/nudus_SRR5767281.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/nudus_SRR5767281.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/pallidus_SRR5767285.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/pallidus_SRR5767285.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/pulcherrimus_DRR107784.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/pulcherrimus_DRR107784.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/pulcherrimus_SRR5767283.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/pulcherrimus_SRR5767283.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/purpuratus_SRR6281818.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/purpuratus_SRR6281818.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/purpuratus_SRR7211988.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/purpuratus_SRR7211988.thresholds.bed.gz'), ('/hb/scratch/mglasena/mosdepth/variegatus_SRR7207203.regions.bed.gz', '/hb/scratch/mglasena/mosdepth/variegatus_SRR7207203.thresholds.bed.gz')]

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
