Questions - What is the difference between species_tree_alignments.py and transcript_create_slt.py?

transcript_create_slt.py was being used to create gene trees for phylonet, so it used the nonparametric bootstrap to make 100 bootstrap trees for each locus. Additionally, it had higher filter thresholds (minimum coverage depth = 15, prop_1x >= 0.9, prop_10x > 0.9, required_gap = 50kb). 

spcies_tree_alignments.py was being used to create gene trees to infer a species tree. No bootstrap trees were needed, and looser filters were applied (minimum coverage depth = 10, prop_1x >= 0.75, prop_10x > 0.75, required_gap = 20kb)