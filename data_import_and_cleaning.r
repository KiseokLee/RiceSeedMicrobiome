####### Sequence data import and cleaning and normalization #######

# import libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metagenomeSeq)
library(vegan)

# import biom file (QIIME2 output) to phyloseq object
bac_phylo=import_biom("OTU_table_final.biom")

# Import sample metadata
map <- read.table(file = 'sample_metadata.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

# Merge biom data object with sample metadata & tree data
phy_tree = read_tree("rooted_tree.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

# Changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", 
                              "Order", "Family", "Genus", "Species")

# Labeling sample names
lab <- c('O.gnu', 'O.bra', 'O.lgl', 'O.rid', 'O.aus', 'O.grd', 'O.alt', 'O.lat', 'O.off', 'O.rhi', 'O.pun', 'O.glu', 'O.lst', 'O.mer', 'O.bar', 'O.gla', 'O.ruf', 'O.niv', 'O.sativa_JS', 'O.sativa_JD', 'O.sativa_IJ', 'O.sativa_PG', 'O.sativa_TJ', 'O.sativa_SS', 'O.sativa_BS', 'O.sativa_CH', 'O.sativa_IP', 'O.sativa_HY', 'O.sativa_SO', 'O.sativa_DC', 'O.sativa_DA', 'O.sativa_DJ', 'O.sativa_DO', 'O.sativa_MY', 'O.sativa_NA', 'O.sativa_PO', 'O.sativa_TI', 'O.sativa_SA', 'O.sativa_GY', 'O.sativa_KN', 'O.sativa_KO', 'O.sativa_AK', 'O.sativa_NP')
sample_data(phy)$Label <- factor(sample_data(phy)$Label , levels =lab)


# Remove host Chloroplast and Mitochondria and Unassigned(Kingdom level) sequences
phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
phy.mt <- subset_taxa(phy, Order == "D_3__Rickettsiales")
phy.un <- subset_taxa(phy, Kingdom == "Unassigned")

# pop_taxa function
pop_taxa = function(physeq, rmTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% rmTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
phy.clean <- pop_taxa(phy, c(vec.cp, vec.mt))

# Remove the "D_3__" patterns for cleaner labels
tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[A-Z]_[0-9]__",replacement="")

# Remove all taxa in negative samples
phy.clean <- pop_taxa(phy.clean, negative_OTUs)

# Remove taxa with sequence length over 300
phy.clean <- pop_taxa(phy.clean, over300 )

# Remove OTUs with less than 5 reads across all samples
v.less5 <- rownames(otu_table(phy.clean)[otu_table(phy.clean) %>% rowSums() < 5])
phy.clean <- pop_taxa(phy.clean, v.less5)

# CSS normalization with metagenomeSeq
met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean)
p <- cumNormStatFast(met.phy.clean)
met.phy.norm <- cumNorm(met.phy.clean, p =p)

# Returns normalized factors for each sample
normFactors(met.phy.norm)
sort(normFactors(met.phy.norm))

# Export normalized count matrix
met.CSS.log <- MRcounts(met.phy.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.phy.norm, norm = T, log = F)


# Convert it back to the phyloseq file
phy.clean.log <- phy.clean
otu_table(phy.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)
phy.clean.nolog <- phy.clean
otu_table(phy.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)