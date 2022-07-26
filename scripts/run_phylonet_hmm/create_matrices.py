import os

# File containing paths to scaffold alignments produced by vcf2phylip. Can create this file using find "$(pwd)" -name "*.nexus" -type f
scaffold_alignments = "scaffold_alignments.txt"
nexus_dir = "/hb/groups/pogson_group/vcf2phylip_4way/prepare_hmm_100_runs/hmm_nexus/"

def add_lines(scaffold_file):
    scaffold_file = str(scaffold_file)
    with open(scaffold_file,'r') as f2:
        inputs = f2.read().splitlines()
        seq = inputs[6]
        length_int = len(seq.split()[1])
        length = str(length_int)
        outdir = "\"" + scaffold_file.split("/")[6] + "\""
        line1 = "#NEXUS"
        line2 = "BEGIN NETWORKS;"
        line3 = "Network net = (intermedius:0.2358047382668817,((fragilis:0.947834903226071,(droebachiensis:0.20102859349321878)#H1:1.0144626112368627::0.578093455107686):0.6877223535737972,(#H1:1.859606432191539::0.42190654489231405,pallidus:1.0045273614964096):0.11154011939943985):4.676738730180861);"
        line4 = "END;"
        line5 = "Begin DATA;"
        line6 = "dimensions ntax=4 nchar=" + length + ";"
        line7 = 'format datatype=dna symbols="ACTG" missing=? gap=-;'
        line8 = "matrix"
        line9 = ";"
        line10 = "END;"
        line11 = "BEGIN PHYLONET;"
        line12 = 'HmmCommand net -gtr -allelemap <intermedius:QB3KMK012; pallidus:QB3KMK002; droebachiensis:QB3KMK014; fragilis:QB3KMK013> -outputdirectory ' + outdir + ' -numberofruns 100 -iterations 1000 -noplots;'
        line13 = "END;"
        with open(scaffold_file.split("/")[6] + ".nexus",'w') as f2: 
            f2.write(line1)
            f2.write("\n")
            f2.write(line2)
            f2.write("\n")
            f2.write(line3)
            f2.write("\n")
            f2.write(line4)
            f2.write("\n")
            f2.write(line5)
            f2.write("\n")
            f2.write(line6)
            f2.write("\n")
            f2.write("\t")
            f2.write(line7)
            f2.write("\n")
            f2.write("\t")
            f2.write(line8)
            f2.write("\n")
            f2.write(inputs[6])
            f2.write("\n")
            f2.write(inputs[7])
            f2.write("\n")
            f2.write(inputs[8])
            f2.write("\n")
            f2.write(inputs[9])
            f2.write("\n")
            f2.write(line9)
            f2.write("\n")
            f2.write(line10)
            f2.write("\n")
            f2.write(line11)
            f2.write("\n")
            f2.write(line12)
            f2.write("\n")
            f2.write(line13)
            f2.write("\n")

def main():
        f = open(scaffold_alignments, "r")
        inputs = f.read().splitlines()
        os.mkdir(nexus_dir)
        os.chdir(nexus_dir)
        for scaffold_file in inputs:
                add_lines(scaffold_file)

if __name__ == "__main__":
        main()
