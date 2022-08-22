import os
from joblib import Parallel, delayed

reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

abba_baba_dir_list = [
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/lividus/all",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/purpuratus/all",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/fragilis/all",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/franciscanus/all",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/nudus/all",
"/hb/home/mglasena/dissertation/data/angsd_abba_baba/depressus/all"
]

def run_abba_baba(dir):
	os.chdir(dir)

	angsd_command = "angsd -doAbbababa 1 -doCounts 1 -baq 1 -ref {} -useLast 1 -bam bam.filelist -out out -blockSize 1000000 -minMapQ 20 -minQ 20 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1".format(reference_genome)

	jackknife = "Rscript jackKnife.R file=out.abbababa indNames=bam.filelistnames outfile=results"

	os.system(angsd_command)
	os.system(jackknife)

def main():
	run_abba_baba()


	Parallel(n_jobs=len(abba_baba_dir_list))(delayed(run_abba_baba)(bam_file) for dir in abba_baba_dir_list)

if __name__ == "__main__":
	main()