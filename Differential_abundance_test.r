######## Differential abundance test with fitZig function from metagenomeSeq #########

## CSS transformation
obj <- phyloseq_to_metagenomeSeq(phy.clean)
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 10, verbose = TRUE)
Type  <-  pData(obj)$Type
mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)

## Differential abundance test with fitZig
ft <- fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
zigFit <- ft$fit
finalMod <- ft$fit$design
contrast.matrix <- makeContrasts(Domesticated - Wild, levels = finalMod)
fit.c <- contrasts.fit(zigFit, contrast.matrix)
fit.e <- eBayes(fit.c)
fit.e
# Get significantly differential taxa
res <- topTable(fit.e, coef=1, adjust="fdr", n=Inf)
resSig <- res[!is.na(res$adj.P.Val), ]
resSig <- data.frame(resSig)

## Construct volcano plot
h = -log10(max(resSig[resSig$adj.P.Val<0.01,]$adj.P.Val))
# FDR Q<0.05 & logFC>|2|, log2FC<|2|
mydata <- resSig %>% mutate(volcano_y = -log10(adj.P.Val))
mydata <- mydata %>% mutate(threshold = ifelse(logFC >= 2 & volcano_y >= h, 'Domesticated', ifelse(logFC <= -2 & volcano_y >= h,"Wild", ifelse(volcano_y < h , "Bottom", 'Middle'))))

theme_set(theme_bw())
ggplot(mydata, aes(x=logFC, y=volcano_y)) +
  xlab('\n Fold Change (Log2)')+
  ylab("-log10 adjusted P\n") +
  geom_point(aes(color = threshold, size=Abun),  alpha=0.6) +
  scale_color_manual(labels = c('Wild'="Wild enriched\n(FDR Q<0.01 & log2FC>2)     ",'Domesticated'= "Domesticated enriched\n(FDR Q<0.01 & log2FC>2)     ",'Middle'='log2FC<|2|     ','Bottom'='NS\n(FDR Q>0.01)     '), values = c("Wild"= "red3", 'Domesticated'='forestgreen',"Middle"="lightblue","Bottom"= "yellow3"))+
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=h, color="maroon4")+
  geom_vline(xintercept=2, color="maroon4")+
  geom_vline(xintercept=0, color="black")+
  geom_vline(xintercept=-2, color="maroon4")+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(-20,20,1))+
  scale_y_continuous(breaks=seq(0,60,10))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
