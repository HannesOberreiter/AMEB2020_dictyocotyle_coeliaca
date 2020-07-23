install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("Biostrings")

library(tidyverse)
library(treeio)
library(ape)
library(ggtree)
library(Biostrings)

btrees <- read.tree("./data/RAxML_bipartitionsBranchLabels.final")
btrees <- ape::root(btrees, c("Schmidtea_mediterranea", "Danio_rerio", "Rattus_norvegicus", "Drosophila_melanogaster", "Cionia_intestinalis"), resolve.root = TRUE)
y <- as.tibble(btrees)

btrees2 <- read.tree("./data/neu3")
btrees2 <- ape::root(btrees2, c("Schmidtea_mediterranea", "Danio_rerio", "Rattus_norvegicus", "Drosophila_melanogaster", "Cionia_intestinalis"), resolve.root = TRUE)


taxa <- list(
  Outgroup = c("Schmidtea_mediterranea", "Danio_rerio", "Rattus_norvegicus", "Drosophila_melanogaster", "Cionia_intestinalis"),
  Trematoda = c("Schistosoma_mansoni", "Clonorchis_sinensis", "Opisthorchis_viverrini"),
  Cestoda = c("Echinococcus_granulosus", "Hymenolepis_microstoma", "Taenia_solium"),
  Monopisthocotylea = c("Leptocotyle_minor", "Dictyocotyle_coeliaca", "Kapentagyrus_tanganicanus", "Gyrodactylus_salaris"),
  Polyopisthocotylea = c("Eudiplozoon_nipponicum", "Diclidophora_denticulata", "Protopolystoma_xenopodis")
)


genetree <- ggtree(btrees2, branch.length="none") + facet_wrap(~.id, ncol=2) + geom_tiplab(size=3) + xlim(0,20)

bootstraptree <- ggtree(x, branch.length = "none") + geom_tiplab(size=3) + xlim(0,20)
bootstraptree

btrees <-groupOTU(btrees, taxa, "Group")
p <- ggtree(btrees, aes(color=Species), branch.length = 'none') + 
  xlim(0, 10) + 
  geom_tiplab(size=3) +
  theme_void() + 
  theme(legend.position="right")
msaplot(p, "data/seq.fas", window=c(50, 500))


x <- as.phylo(x)

groupOTU(bootstraptree, taxa, 'Group') + aes(color=Species) +
  theme(legend.position="right")




test <- read.tree(system.file("extdata/RAxML", "RAxML_bootstrap.H3", package="treeio"))

ggdensitree(rooted, alpha=.3, colour='steelblue', layout = 'fan')+
geom_tiplab(size=3) + xlim(0, 20)

ggdensitree(btrees, alpha=.3, colour='steelblue', branch.length="none")+
  geom_tiplab(size=3) + xlim(0, 20)

s <- c("Clonorchis_sinensis", 
       "Schistosoma_mansoni", 
       "Diclidophora_denticulata", 
       "Eudiplozoon_nipponicum", 
       "Protopolystoma_xenopodis", 
       "Hymenolepis_microstoma", 
       "Taenia_solium",
       "Echinococcus_granulosus",
       "Gyrodactylus_salaris",
       "Kapentagyrus_tanganicanus",
       "Dictyocotyle_coeliaca",
       "Leptocotyle_minor",
       "Schmidtea_mediterranea",
       "Drosophila_melanogaster",
       "Cionia_intestinalis",
       "Rattus_norvegicus",
       "Danio_rerio",
       "Opisthorchis_viverrini"
       )

x <- paste(s, collapse="' | grep -h '")
cat(x)


for(i in s){
  print(i)
  x <- sprintf("grep -l '%s' ./data/genetrees/* | wc -l", i)
  system(x)
}

grep -h 'Clonorchis_sinensis' ./genetrees/* | grep -h 'Schistosoma_mansoni' | grep -h 'Diclidophora_denticulata' | grep -h 'Eudiplozoon_nipponicum' | grep -h 'Protopolystoma_xenopodis' | grep -h 'Hymenolepis_microstoma' | grep -h 'Taenia_solium' | grep -h 'Echinococcus_granulosus' | grep -h 'Gyrodactylus_salaris' | grep -h 'Kapentagyrus_tanganicanus' | grep -h 'Dictyocotyle_coeliaca' | grep -h 'Leptocotyle_minor' | grep -h 'Schmidtea_mediterranea' | grep -h 'Drosophila_melanogaster' | grep -h 'Cionia_intestinalis' | grep -h 'Rattus_norvegicus' | grep -h 'Danio_rerio' | grep -h 'Opisthorchis_viverrini' >> neu3


java -jar DensiTree.jar