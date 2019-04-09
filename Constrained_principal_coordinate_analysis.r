###### Constrained principal coordinate analysis #####

# Import libraries
library(vegan)
library(microbiome)

# Set up dataset
otu <- abundances(phy.clean.log) 
meta <- meta(phy.clean.log)

# CAP analysis #####
# (1) Domestication factor
b.cap.type <- ordinate(phy.clean.log, "CAP", "bray", ~Type)
# Permutest
b.cap.type.perm <- permutest(b.cap.type, permutations = 99999)

# (2) Rice Genome type
phy.genome <- subset_samples(phy.clean.log, Type =='Wild' | Label =='O.sativa_NP' |Label =='O.sativa_TI')
sample_data(phy.genome)$Type <- factor(sample_data(phy.genome)$Type,levels=c('Wild','Domesticated'))
b.cap.genome <- ordinate(phy.genome, "CAP", "bray", ~Genome_type)

# Permutest
b.cap.genome.perm <- permutest(b.cap.genome, permutations = 99999)

# (3) Diversification (Cultivar lineage group of the domesticated rice)
dom_lab <- c('O.sativa_JS', 'O.sativa_JD', 'O.sativa_IJ', 'O.sativa_PG', 'O.sativa_TJ', 'O.sativa_SS', 'O.sativa_BS', 'O.sativa_CH', 'O.sativa_IP', 'O.sativa_HY', 'O.sativa_SO', 'O.sativa_DC', 'O.sativa_DA', 'O.sativa_DJ', 'O.sativa_DO', 'O.sativa_MY', 'O.sativa_NA', 'O.sativa_PO', 'O.sativa_TI', 'O.sativa_SA', 'O.sativa_GY', 'O.sativa_KN', 'O.sativa_KO', 'O.sativa_AK', 'O.sativa_NP')
phy.diversity <- subset_samples(phy.clean.log, Label %in% dom_lab)
b.cap.diversity <- ordinate(phy.diversity, "CAP", "bray", ~Lineage)
# Permutest
b.cap.diversity <- permutest(b.cap.diversity, permutations = 99999)


## Plot CAP
# (1) Domestication factor
plot.b.cap.type     <- plot_ordination(phy.clean.log, b.cap.type, type = "samples",
                                          color = "Type",
                                          shape = "Type")
plot.b.cap.type

plot2.b.cap.type <- plot.b.cap.type + theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+
  xlab('\n Constrained PCoA1 (6.1%)')+
  ylab("MDS1 (11%) \n") +
  ggtitle("Domestication : 6.9% of Variance, P < 0.001 \n") +
  scale_color_manual(values=c('forestgreen','red3'))+
  ## adjust positions
  geom_hline(yintercept=0, color="maroon4",linetype='dotted')+
  geom_vline(xintercept=0, color="maroon4",linetype='dotted')+
  guides(fill = guide_legend(ncol = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(face='bold',size=12, color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot2.b.cap.type

# (2) Genome type - evolution
plot.b.cap.genome    <- plot_ordination(phy.genome, b.cap.genome, type = "samples",
                                       color = "Genome_type",
                                     shape = 'Type')
plot.b.cap.genome

genome_color <- c('firebrick1','darkgreen','dodgerblue3','gold2','lawngreen','cyan1', 'chocolate1'  ,'darkorchid2')

plot2.b.cap.genome <- plot.b.cap.genome + theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+
  xlab('\n Constrained PCoA1 (7.2%)')+
  ylab("Constrained PCoA2 (5.5%) \n") +
  ggtitle("Rice genome type : 24.2% of Variance, P < 0.001 \n") +
  scale_color_manual(values=genome_color)+
  ## adjust positions
  geom_hline(yintercept=0, color="maroon4",linetype='dotted')+
  geom_vline(xintercept=0, color="maroon4",linetype='dotted')+
  guides(fill = guide_legend(ncol = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(face='bold',size=12, color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot2.b.cap.genome

# (3) Diversification (Cultivar lineage group of the domesticated rice)
plot.b.cap.diversity    <- plot_ordination(phy.diversity, b.cap.diversity, type = "samples",
                                        color = "Lineage")
plot.b.cap.diversity

lin_color <- c('firebrick1','darkgreen','dodgerblue3','gold2','lawngreen','cyan1', 'chocolate1'  ,'darkorchid2')

plot2.b.cap.diversity <- plot.b.cap.diversity + theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+
  xlab('\n Constrained PCoA1 (10.3%)')+
  ylab("Constrained PCoA2 (5%) \n") +
  ggtitle("Cultivar lineage group : 24.7% of Variance, P < 0.001 \n") +
  ## adjust positions
  scale_color_manual(values=lin_color)+
  # scale_color_aaas()+
  ## adjust positions
  geom_hline(yintercept=0, color="maroon4",linetype='dotted')+
  geom_vline(xintercept=0, color="maroon4",linetype='dotted')+
  guides(fill = guide_legend(ncol = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(face='bold',size=12, color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot2.b.cap.diversity