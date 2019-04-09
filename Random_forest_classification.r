###### Random forest classification and ridge plot ######

# import libraries
library(phyloseq)
library(dplyr)
library(ggplot2)

random_forest <- function(fun.CSS.nolog, split_ratio, n_tree){
  table <- t(otu_table(fun.CSS.nolog))
  rownames(table)
  k <- 129-75
  result <- c(rep(0,75),rep(1,k))
  table <- cbind(table, result)
  head(table,20)
  dim(table)
  head(table)
  f.dataset <- table
  f.dataset <- data.frame(f.dataset)
  # Encoding the target feature as factor
  f.dataset$result = factor(f.dataset$result, levels = c(0, 1))
  f.dataset$result
  # Splitting the f.dataset into the Training set and Test set
  # install.packages('caTools')
  library(caTools)
  set.seed(123, "L'Ecuyer")
  split = sample.split(f.dataset$result, SplitRatio = split_ratio)
  f.training_set = subset(f.dataset, split == TRUE)
  f.test_set = subset(f.dataset, split == FALSE)
  # Feature Scaling
  ## scaling done with CSS already
  
  # Fitting Random Forest Classification to the Training set
  library(randomForest)
  set.seed(123, "L'Ecuyer")
  f.classifier = randomForest(x = f.training_set[-ncol(f.dataset)],
                              y = f.training_set$result,
                              ntree = n_tree)
  f.training_set[-ncol(f.dataset)]
  f.training_set$result
  
  # Predicting the Test set results
  rf_pred = predict(f.classifier, newdata = f.test_set[-ncol(f.dataset)])
  
  # Making the Confusion Matrix
  cm = table(f.test_set[, ncol(f.dataset)], rf_pred)
  print(cm)
  print(f.classifier)
  # ROC curve
  library(ROCR)
  predictions.rf=as.vector(as.numeric(rf_pred))
  pred.rf=prediction(predictions.rf, f.test_set[,ncol(f.dataset)])
  
  AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
  AUC.rf=AUC.rf@y.values[[1]]
  ROC.rf=performance(pred.rf,"tpr","fpr") #plot the actual ROC curve
  plot(ROC.rf, main="ROC curve", col='green')
  abline(0, 1, col = "grey")
  legend("bottomright", paste0('RandomForest AUC : ',round(as.numeric(AUC.rf), digits = 3)), col = c("green"), pch = c(3))
  
  colnames(f.classifier$err.rate) <- c('OOB','Domestic','Wild')
  plot(f.classifier, main = 'Out-of-bag Error estimate')
  legend("topright", colnames(f.classifier$err.rate),col=1:3,cex=0.8,fill=1:3)
  oob <- f.classifier$err.rate[,1]
  oob.oob <- oob[length(oob)]
  legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob), digits = 4)*100,'%'))
}

random_forest(phy.CSS.log, 0.66, 1500)

## Construct ridge plot

# Create classifier object
table <- otu_table(phy.CSS.log)
table <- t(table)
k <- 129-75
result <- c(rep(0,75),rep(1,k))
table <- cbind(table, result)
dataset <- as.data.frame(table)
dataset$result = factor(dataset$result, levels = c(0, 1))

# Splitting the dataset into the Training set and Test set
library(caTools)
set.seed(123, "L'Ecuyer")
split = sample.split(dataset$result, SplitRatio = 0.66)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

library(randomForest)
set.seed(123, "L'Ecuyer")
classifier = randomForest(x = training_set[-ncol(dataset)],
                          y = training_set$result, ntree=1500)

# Get top 20 discriminant OTU information (by DecreaseMeanGini)
get_imp_tax <- function(phy.clean, classifier){
  imp <- data.frame(importance(classifier))
  imp <- tibble::rownames_to_column(imp, 'OTU')
  imp.desc <- imp %>% arrange(desc(MeanDecreaseGini))
  tax <- tax_table(phy.clean)
  tax <- as.data.frame(tax)
  tax <- tibble::rownames_to_column(tax,'OTU')
  df.clean.ss.5 <- psmelt(phy.clean)
  abun <- df.clean.ss.5 %>% group_by(OTU) %>% summarize(total=sum(Abundance)) %>% arrange(desc(total))
  imp.tax <- left_join(imp.desc, tax, by=c('OTU', 'OTU')) %>% left_join(abun, by=c('OTU','OTU'))
  return(imp.tax)
}

imp_tax <- get_imp_tax(phy.clean, classifier)



# Get dataframe for ridge plot
df.otu <- phy.clean %>% psmelt()
df.otu.rel <- df.otu %>%  
  group_by(Sample) %>%                              
  mutate(RelAbundance = Abundance*1000/sum(Abundance)) 
  
imp.top20 <- imp_tax %>% arrange(desc(MeanDecreaseGini)) %>% head(20)
imp20 <- imp.top20$OTU
df.selected.rel <- df.otu.rel %>% filter(OTU %in% imp20)
df.ridge <- df.selected.rel %>% select(OTU, Sample, RelAbundance)

# Plot ridge plot
plot_ridge <- function(df_ridge){
  ggplot(df_ridge, aes(x=Sample, y=id, height=RelAbundance)) + 
    geom_density_ridges2(stat = "identity", scale=2.5, color='mintcream',size=.3)+
    xlab('\n Rice Species')+
    ylab("Relative Abundance ‰ \n") +
    theme(legend.text=element_text(size=12)) + 
    theme(legend.position="top", legend.spacing.x = unit(0.2, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    geom_vline(xintercept=17.5, color="slategray3", linetype='dashed',size=1)
}

plot_ridge(df.ridge)
