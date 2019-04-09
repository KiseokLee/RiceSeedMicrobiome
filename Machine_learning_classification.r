
###### Classification machine learning testing #######
# import libraries
library(phyloseq)
library(dplyr)
library(ggplot2)

# Data preparation
table <- otu_table(phy.clean.log)
table <- t(table)
k <- 129-75
result <- c(rep(0,75),rep(1,k))
table <- cbind(table, result)
dataset <- as.data.frame(table)

# Encoding the target feature as factor
dataset$result = factor(dataset$result, levels = c(0, 1))
dataset$result

# Splitting the dataset into the Training set and Test set
library(caTools)

set.seed(123)
split = sample.split(dataset$result, SplitRatio = 0.66)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

# (1) Random forest
library(randomForest)
set.seed(123)
classifier = randomForest(x = training_set[-ncol(dataset)],
                            y = training_set$result, ntree=1500)
rf_pred = predict(classifier, newdata = test_set[-ncol(dataset)])
cm = table(test_set[, ncol(dataset)], rf_pred)
print(cm)

## (2) SVM
library(e1071)
classifier_svm = svm(formula = result ~.,
                       data = training_set,
                       type = 'C-classification',
                       kernel ='radial')
svm_pred = predict(classifier_svm, newdata = test_set[-ncol(dataset)])
cm = table(test_set[, ncol(dataset)], svm_pred)
print(cm)

# (3) Naive bayes
library(e1071)
packageVersion('e1071') ## ‘1.7.0’
classifier_nb = naiveBayes(x = training_set[-ncol(dataset)],
                             y = training_set$result)

nb_pred = predict(classifier_nb, newdata = test_set[-ncol(dataset)])
cm = table(test_set[, ncol(dataset)], nb_pred)
print(cm)


# (4) K nearest neighbors
library(class)
packageVersion('class')  # ‘7.3.14’
knn_pred = knn(train = training_set[,-ncol(dataset)], 
               test = test_set[-ncol(dataset)],
               cl = training_set[, ncol(dataset)],
               k=8)
cm = table(test_set[,ncol(dataset)], knn_pred)
cm

# (5) logistic regression
classifier_lr <-  glm(formula = result ~.,
                        family=binomial,
                        data = training_set)
lr_pred <- predict(classifier_lr, type = 'response', newdata=test_set[-ncol(dataset)])
lr_pred <- ifelse(lr_pred > 0.5, 1, 0)
cm = table(test_set[, ncol(dataset)], lr_pred)
print(cm)


### Plotting and comparing ROC curves
library(ROCR)

predictions.rf=as.vector(as.numeric(rf_pred))
pred.rf=prediction(predictions.rf, test_set[,ncol(dataset)])

AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]])
ROC.rf=performance(pred.rf,"tpr","fpr") #plot the actual ROC curve
plot(ROC.rf, main="ROC curve", col='forestgreen')
abline(0, 1, col = "gray")

# Add other ROC curves
roc_addition = function(y_pred, test_set_vector, color, what_method) {
  predictions=as.vector(as.numeric(y_pred))
  pred=prediction(predictions, test_set[,ncol(dataset)])
  AUC=performance(pred,"auc") # Calculate the AUC value
  AUC=AUC@y.values[[1]]
  ROC=performance(pred,"tpr","fpr") # Plot the actual ROC curve
  plot(ROC, main="ROC plot", add=T, col=color)
  my_list <- list('model' = paste0(what_method,': ',round(as.numeric(AUC), digits = 3)), 'color' = color)
  return (my_list)
}

### plot new roc curve easily
roc_addition(svm_pred, test_set[,ncol(dataset)], 'red', 'SVM')
svm <- roc_addition(svm_pred, test_set[,3], 'red', 'SVM') ## for the legend

roc_addition(nb_pred, test_set[,ncol(dataset)], 'blue', 'NaiveBayes')
nb <- roc_addition(nb_pred, test_set[,ncol(dataset)], 'blue', 'NaiveBayes')

roc_addition(knn_pred, test_set[,ncol(dataset)], 'purple', 'K-Nearest neighbors')
knn <- roc_addition(knn_pred, test_set[,ncol(dataset)], 'purple', 'K-Nearest neighbors')

roc_addition(lr_pred, test_set[,ncol(dataset)], 'orange', 'logistic regression')
lr <- roc_addition(lr_pred, test_set[,ncol(dataset)], 'orange', 'logistic regression')


# Add a legend
legend('bottomright',pch =c(3), legend=c(paste0('RandomForest: ',round(as.numeric(AUC.rf), digits = 3)),
                                         svm$model, knn$model,lr$model, nb$model ),
       col=c('forestgreen', svm$color, knn$color,lr$color,  nb$color))




# k-fold validation <- test overfitting 
# (1) random forest cross validation
library(caret)
folds = createFolds(training_set$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set[-x, ]
  test_fold = training_set[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset)],
                              y = training_fold$result, ntree=1500)
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset)])
  cm = table(test_fold[, ncol(dataset)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy
sd <- sd(as.numeric(cv))

# (2) svm cross validation
library(caret)  
folds = createFolds(training_set$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set[-x, ]
  test_fold = training_set[x, ]
  classifier = svm(formula = result ~ .,
                     data = training_fold,
                     type = 'C-classification',
                     kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset)])
  cm = table(test_fold[, ncol(dataset)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ## 0.8972222
sd <- sd(as.numeric(cv))

# (3) nb cross validation
library(caret)
folds = createFolds(training_set$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set[-x, ]
  test_fold = training_set[x, ]
  classifier = naiveBayes(x = training_fold[-ncol(dataset)],
                            y = training_fold$result)
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset)])
  cm = table(test_fold[, ncol(dataset)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   
sd(as.numeric(cv)) 

# (4) knn cross validation
library(caret)
folds = createFolds(training_set$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set[-x, ]
  test_fold = training_set[x, ]
  knn_pred = knn(train = training_fold[,-ncol(dataset)], 
                 test = test_fold[-ncol(dataset)],
                 cl = training_fold[, ncol(dataset)],
                 k=8)
  # Making the Confusion Matrix
  cm = table(test_fold[,ncol(dataset)], knn_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   
sd(as.numeric(cv)) 

# (5) logistic regression cross validation
library(caret)
folds = createFolds(training_set$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set[-x, ]
  test_fold = training_set[x, ]
  classifier_lr <-  glm(formula = result ~.,
                          family=binomial,
                          data = training_fold)
  lr_pred <- predict(classifier_lr, type = 'response', newdata=test_fold[-ncol(dataset)])
  lr_pred <- ifelse(lr_pred > 0.5, 1, 0)
  # Making the Confusion Matrix
  cm = table(test_fold[, ncol(dataset)], lr_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   
sd(as.numeric(cv)) 
