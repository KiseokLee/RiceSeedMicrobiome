
## (1) Phylosymbiosis analysis for Bacterial Kingdom
##Dendrogram log transformation
SeedSp.m <- otu_table(phy.clean.log.m) #merged normalized first part
SeedSp.m

#mat = t(otu.43g) #for non merged OTU table

mat = t(SeedSp.m)
d = (vegdist(mat, method="bray"))
h = hclust(d)

plot(h, main = "Bray Curtis method", sub = "", xlab="", axes = FALSE, hang = -1)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))

tree.ours.m <- phy_tree(phy.clean.log.m)
plot(tree.ours.m)

host_tree <- read.tree("Newick Export.nwk")
host_tree$tip.label

plot(host_tree)
RFmeasures(host_tree,h,nRandom=100000)

Statistics<-Phylosymbiosis(t(SeedSp.m),host_tree,tree.ours.m,gen.mat, TraitDist=NA)
Statistics

###Extract tree file
#Compute Beta
unifracs <- aperm(GUniFrac(t(SeedSp.m), tree.ours.m, alpha=c(0, 0.5, 1))$unifracs,c(3,1,2)) #Unifrac
bray <-vegdist(x=t(SeedSp.m),method="bray") #vegan
jac <-vegdist(x=t(SeedSp.m),method="jaccard") #vegan


#Compute correlation
# Original method
# Hierchachical clustering of microbiomes
Dendro_W_Uni=hclust(as.dist(unifracs['d_1',,]),method = "average")
Dendro_UW_Uni=hclust(as.dist(unifracs['d_UW',,]),method = "average")
Dendro_Bray=hclust(bray,method = "average")
Dendro_Jac=hclust(jac,method = "average")

class(Dendro_W_Uni) # must be hclust class
my_tree_wuni <- as.phylo(Dendro_W_Uni) 
write.tree(phy=my_tree_wuni, file="Dendro_W_Uni.newick")

class(Dendro_UW_Uni) # must be hclust class
my_tree_uw_uni <- as.phylo(Dendro_UW_Uni) 
write.tree(phy=my_tree_uw_uni, file="Dendro_UW_Uni.newick")

class(Dendro_Bray) # must be hclust class
my_tree_bray <- as.phylo(Dendro_Bray) 
write.tree(phy=my_tree_bray, file="Dendro_Bray.newick")

class(Dendro_Jac) # must be hclust class
my_tree_jac <- as.phylo(Dendro_Jac) 
write.tree(phy=my_tree_jac, file="Dendro_Jac.newick")


###Cophylogenetic analysis (PACO)
install.packages('paco')
library(paco)
library(ape)

rice_tree<-cophenetic(host_tree)
bac_tree<- as.matrix(cophenetic(Dendro_Bray))
bac_tree
association <- cbind(host_tree$tip.label, host_tree$tip.label)
association
cophyloplot(host_tree, my_tree, assoc = association,
            length.line = 4, space = 28,gap = 3, rotate = FALSE,
            col = par("fg"),
            show.tip.label = TRUE, font = 3)

D_bac <-prepare_paco_data(H=rice_tree, P=bac_tree, HP = gl_links)

## (2) Phylosymbiosis analysis for Fungal Kingdom
##Dendrogram log transformation
SeedSp.m <- otu_table(fun.clean.log.m) #merged normalized first part
SeedSp.m

#mat = t(otu.43g) #for non merged OTU table

mat = t(SeedSp.m)
d = (vegdist(mat, method="bray"))
h = hclust(d)
plot(h, main = "Bray Curtis method", sub = "", xlab="", axes = FALSE, hang = -1)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))

tree.ours.m <- phy_tree(fun.clean.log.m)
plot(tree.ours.m)

host_tree <- read.tree("Newick Export.nwk")
host_tree$tip.label

plot(host_tree)
RFmeasures(host_tree,h,nRandom=100000)

genetic.distance <- read.table('Oryza_genus_wTI_distance_matrix.tsv')
gen.mat <- as.matrix(genetic.distance)

Statistics<-Phylosymbiosis(t(SeedSp.m),host_tree,tree.ours.m,gen.mat, TraitDist=NA)
Statistics

###Extract tree file
#Compute Beta
unifracs <- aperm(GUniFrac(t(SeedSp.m), tree.ours.m, alpha=c(0, 0.5, 1))$unifracs,c(3,1,2)) #Unifrac
bray <-vegdist(x=t(SeedSp.m),method="bray") #vegan
jac <-vegdist(x=t(SeedSp.m),method="jaccard") #vegan


#Compute correlation
# Original method
# Hierchachical clustering of microbiomes
Dendro_W_Uni=hclust(as.dist(unifracs['d_1',,]),method = "average")
Dendro_UW_Uni=hclust(as.dist(unifracs['d_UW',,]),method = "average")
Dendro_Bray=hclust(bray,method = "average")
Dendro_Jac=hclust(jac,method = "average")


class(Dendro_W_Uni) # must be hclust class
my_tree_wuni <- as.phylo(Dendro_W_Uni) 
write.tree(phy=my_tree_wuni, file="Dendro_W_Uni.newick")

class(Dendro_UW_Uni) # must be hclust class
my_tree_uw_uni <- as.phylo(Dendro_UW_Uni) 
write.tree(phy=my_tree_uw_uni, file="Dendro_UW_Uni.newick")

class(Dendro_Bray) # must be hclust class
my_tree_bray <- as.phylo(Dendro_Bray) 
write.tree(phy=my_tree_bray, file="Dendro_Bray.newick")

class(Dendro_Jac) # must be hclust class
my_tree_jac <- as.phylo(Dendro_Jac) 
write.tree(phy=my_tree_jac, file="Dendro_Jac.newick")





