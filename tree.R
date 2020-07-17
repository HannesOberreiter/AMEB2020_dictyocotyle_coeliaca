#install.packages("TreeTools")
#install.packages("apex")
library(apex)
library(TreeTools)
library(here)

NEWICK <- ape::read.tree("./data/RAxML_bipartitionsBranchLabels.final")
NEWICK <- root(NEWICK, outgroup = c("Schmidtea_mediterranea"))
NEWICK$Nnode
plot.phylo(NEWICK, type="phylogram",  use.edge.length = FALSE)

NEWICKdev.off()