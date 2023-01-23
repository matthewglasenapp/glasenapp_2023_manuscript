import os
import multiprocessing
from joblib import Parallel, delayed

num_scaffolds = 21

root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/"

# Directory for phylonet_hmm scaffold input files 
output_dir = root_dir + "/hmm_nexus_files/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

# Copy/paste from phylonet InferNetwork_ML run with 1 retiulcation age
phylogenetic_network = "Network net = (pulcherrimus,((#H1,pallidus),(fragilis,(droebachiensis)#H1)));"
allele_map = "<pulcherrimus:QB3KMK016; pallidus:QB3KMK002; droebachiensis:QB3KMK014; fragilis:QB3KMK013>"

number_taxa = 4 
number_runs = 100
number_iterations = 1000
threads = 6

def get_scaffold_file_paths():
    # Create file containing paths to scaffold alignments produced by vcf2phylip.
    create_scaffold_alignment_paths_file = 'find {} "$(pwd)" -name "*.nexus" -type f > scaffold_alignment_paths_file'.format(root_dir)
    os.system(create_scaffold_alignment_paths_file)

    with open("scaffold_alignment_paths_file", "r") as f:
        scaffold_alignment_file_path_list = f.read().splitlines()

    os.system("rm scaffold_alignment_paths_file")
    return scaffold_alignment_file_path_list

def create_hmm_input_file(scaffold_file):
    with open(str(scaffold_file),'r') as f:
        inputs = f.read().splitlines()
    seq = inputs[6]
    length = str(len(seq.split()[1]))
    scaffold_name = scaffold_file.split("/")[-1].split(".min4")[0]
    outdir = "\"" + scaffold_name + "\""
    line1 = "#NEXUS"
    line2 = "BEGIN NETWORKS;"
    line3 = phylogenetic_network
    line4 = "END;"
    line5 = "Begin DATA;"
    line6 = "dimensions ntax={} nchar={};".format(number_taxa,length)
    line7 = 'format datatype=dna symbols="ACTG" missing=? gap=-;'
    line8 = "matrix"
    line9 = ";"
    line10 = "END;"
    line11 = "BEGIN PHYLONET;"
    line12 = "HmmCommand net -gtr -allelemap {} -outputdirectory {} -numberofruns {} -iterations {} -threads {} -noplots;".format(allele_map, outdir, number_runs, number_iterations, threads)
    line13 = "END;"
    

    with open(output_dir + scaffold_name + ".nexus",'w') as f2: 
        f2.write(line1 + "\n")
        f2.write(line2 + "\n")
        f2.write(line3 + "\n")
        f2.write(line4 + "\n")
        f2.write(line5 + "\n")
        f2.write(line6 + "\n" + "\t")
        f2.write(line7 + "\n" + "\t")
        f2.write(line8 + "\n")
        f2.write(inputs[6] + "\n")
        f2.write(inputs[7] + "\n")
        f2.write(inputs[8] + "\n")
        f2.write(inputs[9] + "\n")
        f2.write(line9 + "\n")
        f2.write(line10 + "\n")
        f2.write(line11 + "\n")
        f2.write(line12 + "\n")
        f2.write(line13 + "\n")

def main():
        scaffold_alignment_file_path_list = get_scaffold_file_paths()

        Parallel(n_jobs=num_scaffolds)(delayed(create_hmm_input_file)(scaffold_alignment_file) for scaffold_alignment_file in scaffold_alignment_file_path_list)

if __name__ == "__main__":
        main()
