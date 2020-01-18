## (1) Phylosymbiosis analysis for Bacterial Kingdom
library(parallel)
library(geiger)
library(abind)
library(vegan)
library(phytools)
library(GUniFrac)
library(ape)
source(PhyloSymbiosis_Functions.R) 
#Available at https://github.com/FloMazel/Phylosymbiosis-Ecological-model/blob/master/PhyloSymbiosis_Functions.R

### Use Wild and 3 domesticated rice
phy.clean.log.first
SeedSp.m <- otu_table(phy.clean.log.first)

#Rice phylogenetic tree
host_tree <- read.tree("Newick Export.nwk")
host_tree$tip.label
plot(host_tree)

#Compute Beta
bray <-vegdist(x=t(SeedSp.m),method="bray") #vegan

#Compute correlation
# Original method
# Hierchachical clustering of microbiomes
Dendro_Bray=hclust(bray,method = "average")

class(Dendro_Bray) # must be hclust class
my_tree_bray <- as.phylo(Dendro_Bray) 
write.tree(phy=my_tree_bray, file="Dendro_Bray.newick")

##Cophylogenetic tree for our data
library(phytools)
my_tree_bray# Dendrogram
host_tree # rice phylogeny

cophylo_bac<-cophylo(host_tree,my_tree_bray, rotate = F)
plot(cophylo_bac)
plot.cophylo(cophylo_bac, scale.bar=c(0.1,0.1))

##Mantel test
genetic.distance <- read.table('Oryza_genus_wTI_distance_matrix.tsv')
gen.mat <- as.matrix(genetic.distance)
mantel.bray<-mantel(bray, gen.mat, method = "spearman")
RFmeasures(host_tree,SymbiontDendro=Dendro_Bray,nRandom=100000)


## (2) Phylosymbiosis analysis for Fungal Kingdom
# Use Wild and 3 domesticated rice
fun.clean.log.first
SeedSp.m <- otu_table(fun.clean.log.first)

#Rice phylogenetic tree
host_tree <- read.tree("Newick Export.nwk")
host_tree$tip.label
plot(host_tree)

#Compute Beta
bray <-vegdist(x=t(SeedSp.m),method="bray") #vegan

#Compute correlation
# Original method
# Hierchachical clustering of microbiomes
Dendro_Bray=hclust(bray,method = "average")

class(Dendro_Bray) # must be hclust class
my_tree_bray <- as.phylo(Dendro_Bray) 
write.tree(phy=my_tree_bray, file="Dendro_Bray.newick")

##Cophylogenetic tree for our data
library(phytools)
my_tree_bray# Dendrogram
host_tree # rice phylogeny

cophylo_bac<-cophylo(host_tree,my_tree_bray, rotate = F)
plot(cophylo_bac)
plot.cophylo(cophylo_bac, scale.bar=c(0.1,0.1))

##Mantel test
genetic.distance <- read.table('Oryza_genus_wTI_distance_matrix.tsv')
gen.mat <- as.matrix(genetic.distance)
mantel.bray<-mantel(bray, gen.mat, method = "spearman")
RFmeasures(host_tree,SymbiontDendro=Dendro_Bray,nRandom=100000)












