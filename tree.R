#install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("treeio")
install.packages("tidyverse")
install.packages("patchwork")

library(tidyverse)
library(patchwork)
library(treeio)
library(ggtree)
library(ape)

# our taxa List
TAXA <- list(
  Outgroup = c("Schmidtea mediterranea", "Danio rerio", "Rattus norvegicus", "Drosophila melanogaster", "Cionia intestinalis"),
  Trematoda = c("Schistosoma mansoni", "Clonorchis sinensis", "Opisthorchis viverrini"),
  Cestoda = c("Echinococcus granulosus", "Hymenolepis microstoma", "Taenia solium"),
  Monop = c("Leptocotyle minor", "Dictyocotyle coeliaca", "Kapentagyrus tanganicanus", "Gyrodactylus salaris"),
  Polyp = c("Eudiplozoon nipponicum", "Diclidophora denticulata", "Protopolystoma xenopodis")
)
CLADES <- names(TAXA)
COLORS <- c("grey", "black", "orange", "darkblue", "blue")

# load tree
FINAL     <- treeio::read.newick("./data/RAxML_bipartitions.final", node.label = "support")	
# rename (remove underscore)
RENAME    <- tibble(org = FINAL@phylo[["tip.label"]], new = FINAL@phylo[["tip.label"]] %>% str_replace_all("_", " "))
FINAL     <- rename_taxa(FINAL, RENAME)
# transfer to phylo for ggtree, TODO we may can use the tree format but don't know how to load bootstrap
FINALTREE <- treeio::as.phylo(FINAL)
FINALTREE[["node.label"]] <- FINAL@data[["support"]]
# root
FINALTREE <- treeio::root(FINALTREE, TAXA$Outgroup)
# create clades
FINALTREE <- groupOTU(FINALTREE, TAXA, group_name = "Clades")
# fetch our clade nodes
CLADENODES <- map(TAXA, function(x){
  MRCA(FINALTREE, x)
})

# Our Tree
MAIN_TREE <- ggtree::ggtree(FINALTREE, branch.length="none") +
  geom_hilight(23, alpha=0.2) +
  geom_hilight(29, alpha=0.2) +
  theme(legend.position = 'none')+
  geom_tiplab(size=3, aes(color=Clades), fontface='italic') + 
  geom_label2(
    aes(
      label=label, 
      subset=(!is.na(as.numeric(label) | label=="Root"))), 
    color = "black", size=2
    )+
  #geom_nodepoint() +
  #geom_nodelab(size=1.75, hjust = 2, vjust = -0.3) + 
  #geom_rootedge(rootedge = 0.1) +
  scale_color_manual(breaks=CLADES, values=COLORS) +
  xlim(NA,10)

# colour lower value bootstrap
d <- MAIN_TREE$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label < 100,]
MAIN_TREE = MAIN_TREE + geom_label(data=d, aes(label=label), color = "red", size=2)

for(i in 1:length(CLADENODES)){
  dnode = CLADENODES[[i]]
  dname = CLADES[i]
  MAIN_TREE = MAIN_TREE +
   geom_cladelabel(
     node=dnode, label=dname, align=T, 
     fontsize=3, angle=90, offset=3,
     offset.text=0.2, hjust="center", color=COLORS[i]
     )
}

## Plot our main Tree
MAIN_TREE
ggplot2::ggsave("./data/tree.pdf", MAIN_TREE)

## get the gene tree

# create original file names
TAXA_ORG <- lapply(TAXA, function(x){str_replace_all(x," ", "_")})
TAXA_LIST <- TAXA_ORG %>% unlist()

# print number busco genes for each taxa
for(i in TAXA_LIST){
  print(i)
  x <- sprintf("grep -l '%s' ./data/genetrees/* | wc -l", i)
  system(x)
}

# only get newick trees were genes are available for all taxa
#x <- paste(TAXA_ORG, collapse="' | grep -h '")
#x <- paste("grep -h -f ./data/genetrees/* -e '", x, "' > ./data/allbusco.tree", sep ="")
#system(x)
# only get newick trees were genes are available for all taxa + filename
x <- paste(TAXA_ORG, collapse="' | grep '")
x <- paste("grep -f ./data/genetrees/* -e '", x, "' > ./data/allbusco_names.tree", sep ="")
system(x)

# read our trees
BUSCOTREES <- read.tree("./data/allbusco_names.tree")
# change names to reflect genes
names(BUSCOTREES) <- names(BUSCOTREES) %>% 
  str_extract("EOG(.*).clustalo") %>% str_replace(".clustalo", "")
GENE_NAMES <- names(BUSCOTREES)

PLOT_LIST = list()

for(b in 1:length(BUSCOTREES)){
  gene <- GENE_NAMES[b]
  print(gene)
  GENE_TREE <- treeio::as.phylo(BUSCOTREES[[b]])
  GENE_TREE <- rename_taxa(GENE_TREE, RENAME)
  GENE_TREE <- treeio::root(GENE_TREE, TAXA$Outgroup)
  GENE_TREE <- groupOTU(GENE_TREE, TAXA, group_name = "Clades")
  CLADENODES <- map(TAXA, function(x){
    MRCA(GENE_TREE, x)
  })
  PLOT_TREE <- ggtree::ggtree(GENE_TREE, branch.length="none") +
    ggtitle(paste("Gene: ", gene))+
    theme(legend.position = 'none')+
    geom_tiplab(size=3, aes(color=Clades), fontface='italic') + 
    scale_color_manual(breaks=CLADES, values=COLORS) +
    xlim(NA,14)
  
  for(i in 1:length(CLADENODES)){
    if(i == 1){next}
    dnode = CLADENODES[[i]]
    dname = CLADES[i]
    PLOT_TREE = PLOT_TREE +
      geom_cladelabel(
        node=dnode, label=dname, align=T, 
        fontsize=3, angle=90, offset=3.5,
        offset.text=0.2, hjust="center", color=COLORS[i]
      )
  }

  PLOT_LIST[[ gene ]] <- PLOT_TREE
}


p <- (PLOT_LIST[[1]] | PLOT_LIST[[2]]) / (PLOT_LIST[[3]] | PLOT_LIST[[4]])
ggplot2::ggsave("./data/gene_tree.pdf", p, units = "cm", width=35, height=20)
 








