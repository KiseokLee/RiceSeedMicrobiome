###### Alpha diversity and Lorenz curve #####

# Import libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(reshape2)
library(vegan)



# get the number of observed otus for each labels
get_obs_b <- function(phy.clean.nolog){
  df.otu <- data.frame(t(otu_table(phy.clean.nolog)))
  df.df.otu <- tibble::rownames_to_column(df.df.otu, var = 'SampleID')
  v <- df.df.otu[1,]
  df.otu$Observed <- apply(df.otu[-1], 1, function(x) sum(x!= 0))
  df.obs <- df.otu[,c(1,length(df.otu))]
  return(df.obs)
}

# get diversity indices dataframe
get_diversity_df <- function(phy.clean.nolog){
  df.obs <- get_obs_b(phy.clean.nolog)
  df.otu <- t(otu_table(phy.clean.nolog))
  shannon.nolog <- vegan::diversity(df.otu, index='shannon')
  invsimpson.nolog <- vegan::diversity(df.otu, index='invsimpson')
  df.shannon.nolog <- data.frame(Shannon=shannon.nolog)
  df.invsimpson.nolog <- data.frame(invSimpson=invsimpson.nolog)
  df.shannon.nolog <- tibble::rownames_to_column(df.shannon.nolog, var = 'SampleID')
  df.invsimpson.nolog  <- tibble::rownames_to_column(df.invsimpson.nolog , var = 'SampleID')
  meta <- data.frame(sample_data(phy.clean.nolog)[,c('Type','Label'),])
  df.meta  <- tibble::rownames_to_column(meta , var = 'SampleID')
  df_diversity.nolog <- df.obs %>% left_join(df.shannon.nolog,by=c('SampleID'='SampleID')) %>% left_join(df.invsimpson.nolog,by=c('SampleID'='SampleID')) %>% left_join(df.meta,by=c('SampleID'='SampleID'))
  return(df_diversity.nolog )
}

df.diversity.nolog <- get_diversity_df(phy.clean.nolog)

# Hellinger transformation diversity dataframe
phy.clean.hell = transform_sample_counts(phy.clean, function(x){sqrt(x / sum(x))})
df.diversity.hell <- get_diversity_df(phy.clean.hell)

## Comparing wild and domesticated rice alpha diversity 
# (1) Number of observed OTUs
par(font.axis = 2, font.lab=2)
boxplot(Observed ~ Type, data = df.diversity.hell, ylab = "Observed OTUs", outline=F, xaxt="n")
# Testing whether it has parametric distribution
histogram(subset(df.diversity.hell, Type=='Domesticated')$Observed)
histogram(subset(df.diversity.hell, Type=='Wild')$Observed)
# Wilcoxon test
x <- subset(df.diversity.hell, Type=='Domesticated')$Observed
y <- subset(df.diversity.hell, Type=='Wild')$Observed
wilcox.test(x, y, conf.int = TRUE) 

# (2) InvSimpson Index
boxplot(invSimpson ~ Type, data = df.diversity.hell, ylab = "Inverse Simpson Index", outline=F, xaxt="n")
# Testing whether it has parametric distribution
histogram(subset(df.diversity.hell, Type=='Domesticated')$invSimpson)
histogram(subset(df.diversity.hell, Type=='Wild')$invSimpson)
# wilcoxon test
x <- subset(df.diversity.hell, Type=='Domesticated')$invSimpson
y <- subset(df.diversity.hell, Type=='Wild')$invSimpson
wilcox.test(x, y, conf.int = TRUE)

# (3) Shannon Index
boxplot(Shannon ~ Type, data = df.diversity.hell, ylab = "Shannon Index", outline=F)
# parametric?
histogram(subset(df.diversity.hell, Type=='Domesticated')$Shannon)
histogram(subset(df.diversity.hell, Type=='Wild')$Shannon)
# wilcoxon test
x <- subset(df.diversity.hell, Type=='Domesticated')$Shannon
y <- subset(df.diversity.hell, Type=='Wild')$Shannon
wilcox.test(x, y, conf.int = TRUE)
