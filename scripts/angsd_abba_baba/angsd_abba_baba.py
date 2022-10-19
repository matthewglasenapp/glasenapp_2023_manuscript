import os

reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

threads = 4

abba_baba_dir_list = [
#"/hb/home/mglasena/dissertation/data/angsd_abba_baba/lividus/all/",
#"/hb/home/mglasena/dissertation/data/angsd_abba_baba/purpuratus/all/",
#"/hb/home/mglasena/dissertation/data/angsd_abba_baba/fragilis/all/",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/franciscanus/all/",
#"/hb/home/mglasena/dissertation/data/angsd_abba_baba/nudus/all/",
#"/hb/home/mglasena/dissertation/data/angsd_abba_baba/depressus/all/"
]

def run_abba_baba(dir):
	os.chdir(dir)

	angsd_command = "angsd -doAbbababa 1 -doCounts 1 -baq 1 -ref {} -useLast 1 -bam bam.filelist -out out_mq30 -blockSize 1000000 -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -nThreads {}".format(reference_genome, threads)

	jackknife = "Rscript jackKnife.R file=out_mq30.abbababa indNames=bam.filelistnames outfile=results_mq30"

	os.system(angsd_command)
	os.system(jackknife)

def main():
	#array_id = os.environ["array_id"]
	#print("Array ID: {}".format(array_id))
	#dir = abba_baba_dir_list[int(array_id)]
	#print("Working in {}".format(dir))

	run_abba_baba("/hb/home/mglasena/dissertation/data/angsd_abba_baba/franciscanus/all/")

	print("Analysis Finished!")

if __name__ == "__main__":
	main()