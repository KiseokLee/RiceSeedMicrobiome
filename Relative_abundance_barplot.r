###### Relative abundance barplot with phyloseq & ggplot ######

# import libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(reshape2)
library(forcats) 
library(magrittr)


# Melt phyloseq object into dataframe
df.genus <- phy %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

# Set horizontal order of barplot
lb <- c('O.gnu', 'O.bra', 'O.lgl', 'O.rid', 'O.aus', 'O.grd', 'O.alt', 'O.lat', 'O.off', 'O.rhi', 'O.pun', 'O.glu', 'O.lst', 'O.mer', 'O.bar', 'O.gla', 'O.ruf', 'O.niv', 'O.sativa_JS', 'O.sativa_JD', 'O.sativa_IJ', 'O.sativa_PG', 'O.sativa_TJ', 'O.sativa_SS', 'O.sativa_BS', 'O.sativa_CH', 'O.sativa_IP', 'O.sativa_HY', 'O.sativa_SO', 'O.sativa_DC', 'O.sativa_DA', 'O.sativa_DJ', 'O.sativa_DO', 'O.sativa_MY', 'O.sativa_NA', 'O.sativa_PO', 'O.sativa_TI', 'O.sativa_SA', 'O.sativa_GY', 'O.sativa_KN', 'O.sativa_KO', 'O.sativa_AK', 'O.sativa_NP')
df.genus$Sample <- factor(df.genus$Sample, levels = lb)

# Make NA to unidentified
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))

# Label taxa with less than 5% relative abundance to 'Low abundance'
levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')
df.genus.rel <- df.genus %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 5,]$Genus <- 'Low abundance'

# Set vertical order of genus by total abunance
gen.ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.gen <- gen.ord$Genus
df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.gen)

# Plot relative abundance
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=Sample, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack', colour="black") +
  scale_fill_manual(values = genus_color) +  
  xlab('')+
  ylab("Relative Abundance ‰ \n") +
  guides(fill = guide_legend(ncol = 7,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank()) 

df.genus.rel.p1