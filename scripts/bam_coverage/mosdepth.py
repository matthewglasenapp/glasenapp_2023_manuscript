regions_file = "regions.bed"
threads = 4

def get_bam_file_paths():
	get_files = 'find /hb/scratch/mglasena/data/dedup_mapped_bam_files/ -type f -name "*.bam*" | grep -v "bai" > bam_file_paths.txt'
	os.system(get_files)

	with open("bam_file_paths.txt",r) as f:
		bam_file_paths_list = f.read().splitlines()

	return bam_file_paths_list

def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by {} --no-per-base -t {} --fast-mode {} {}".format(regions_file, threads, prefix, bam_file)
	os.system(mosdepth)

	global_dist_file = prefix + ".mosdepth.global.dist.txt"
	dist_outfile = prefix + ".dist.html"
	plot_dist = "python3 plot-dist.py --output {} {}".format(dist_outfile, global_dist_file)
	os.system(plot_dist)

def main():
	bam_file_paths_list = get_bam_file_paths()
	for bam_file in bam_file_paths_list:
		run_mosdepth(bam_file)
