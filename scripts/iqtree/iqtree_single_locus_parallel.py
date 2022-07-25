import os
import multiprocessing
from joblib import Parallel, delayed

num_cores = multiprocessing.cpu_count()

fasta_alignment_directory = ""
fasta_alignment_file_list = os.listdir(fasta_alignment_directory)

def iqtree(file):
    run_iqtree = "iqtree -s {} -m MFP".format(file)
    os.system(run_iqtree)

main():
    Parallel(n_jobs=num_cores)(delayed(iqtree)(fasta_alignment_file) for fasta_alignment_file in fasta_alignment_file_list)
    
if __name__ == "__main__":
    main()
    