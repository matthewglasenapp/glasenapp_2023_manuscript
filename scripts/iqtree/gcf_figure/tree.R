library(ggtree)
library(phytools)
library(treeio)
library(svglite)

# Parse IQ-TREE files
#tree_data <- read.iqtree("gcf.txt")
#tree_data

# Species tree with no branch lengths
#tree <- read.tree("no_bl.txt")

# Species tree with branch lengths 
#tree <- read.tree("gcf.txt")

#rooted_tree = midpoint.root(tree)
#ape::write.tree(rooted_tree)

# Species tree rooted with midpoint.root(), then manually add the branch support values
rooted_tree <- read.tree("gcf.treefile")

# Display internal node numbers
#ggtree(rooted_tree) + geom_text(aes(label=node), hjust=-.3)

par(mar = c(5.1, 4.1, 4.1, 2.1))

figure <- ggtree(rooted_tree, ladderize=FALSE) + geom_tiplab(size=6) + geom_cladelab(node=11, label="S Clade", align=T, offset = .018, textcolor='purple', barcolor='purple', barsize = 1, fontsize = 5) + geom_cladelab(node=16, label="M Clade", align=T, offset = 0.013, textcolor='red', barcolor='red', barsize = 1, fontsize = 5) + geom_treescale() + xlim(0,0.05) + geom_text2(aes(subset = !isTip, label=label), hjust=0, vjust=-.5, size=4) 

figure

ggsave(plot = figure, "test_tree.svg", width=169, height = 150, units = "mm")

# Highlighting clades
#+ geom_hilight(node=11, fill="purple", alpha=.4) + geom_hilight(node=16, fill="red", alpha=.4)
