#############################################
# Load the required libraries and functions #
#############################################
library(survminer)
library(survival)
library(randomForestSRC)
library(stringr)
library(dplyr)
library(survival)
fitModel <- function(survivalTime,status,features){
  data <- list(time=survivalTime,status=status,x=features)
  fit <- coxph(Surv(time,status)~x, data)
  return(fit)
}

final_selected_features <- readRDS(file.path(normalizePath("."),"selected_phenotypes.rds"))
#
#for model, will use frequency of CD154- phenotypes in unstimulated samples and
#the difference in frequency of CD154+ phenotypes between stimulated and unstimulated samples
#
final_selected_features <- c(paste0("U_",final_selected_features[grepl("CD154-",final_selected_features)]),
                             paste0("S_",final_selected_features[grepl("CD154\\+",final_selected_features)]))

#
#bring in counts from training set, transform to frequencies
#
countDF <- as.data.frame(readRDS(file.path(normalizePath("."),"part03","faustData",
                                           "faustCountMatrix.rds")))
parentCount <- apply(countDF,1,sum)
allCellPops <- setdiff(colnames(countDF),"0_0_0_0_0")
propDF <- countDF[,allCellPops]
for (pCellPop in allCellPops) {
    obsCounts <- propDF[,pCellPop]
    propDF[,pCellPop] <- ((obsCounts)/(parentCount))
} 
propDF$name <- as.character(rownames(propDF))

#
#merge on meta data
#
mpFixed <- allCellPops
metaData <- readRDS(file.path(normalizePath("."),"flowReMi_metaData.rds"))
allFeatureDF <- inner_join(propDF,metaData,by=c("name"))
allFeatureDF <- allFeatureDF[,-which(colnames(allFeatureDF) %in% c("fileNum","name"))]

#
#subset to stim samples
#
stimDF <- allFeatureDF[which(allFeatureDF$stimStatus==1),]
stimDF <- stimDF[,-which(colnames(stimDF)=="stimStatus")]
colnames(stimDF)[which(colnames(stimDF) %in% allCellPops)] <- paste0("S_",colnames(stimDF)[which(colnames(stimDF) %in% allCellPops)])

#
#subset to unstim samples
#
unstimDF <- allFeatureDF[which(allFeatureDF$stimStatus==0),]
unstimDF <- unstimDF[,-which(colnames(unstimDF)=="stimStatus")]
colnames(unstimDF)[which(colnames(unstimDF) %in% allCellPops)] <- paste0("U_",colnames(unstimDF)[which(colnames(unstimDF) %in% allCellPops)])

#
#combine, derive differences
#
trainDF <- inner_join(stimDF,unstimDF,by=c("subjectID","Status","Survival.Time","inTraining"))
for (cellPop in allCellPops) {
    trainDF[,paste0("S_",cellPop)] <- (trainDF[,paste0("S_",cellPop)]-trainDF[,paste0("U_",cellPop)])
}
trainDF$Survival.Time <- trainDF$Survival.Time-min(trainDF$Survival.Time)+1

#
#fit the model on the training set
#
inputDF <- trainDF[,c(final_selected_features,"Survival.Time","Status"),drop=FALSE]
rf_training_fit <- rfsrc(
    formula = Surv(Survival.Time,Status)~.,
    data = inputDF,
    ntree = 500,
    forest = T,
    na.action = "na.impute",
    seed = 12345
)
to_predict_training <- trainDF[,final_selected_features,drop=FALSE]
prediction_training_set <- predict(rf_training_fit,newdata=to_predict_training, na.action="na.impute")
train_survivalTime <- trainDF[,"Survival.Time"]
train_status <- trainDF[,"Status"]
fitted_values_training_set <- fitModel(train_survivalTime,train_status,prediction_training_set$predicted)
pvalue_training_set <- summary(fitted_values_training_set)$logtest[3]
cindex_training_set <- summary(fitted_values_training_set)$concordance
print(paste0("-log10(pvalue) on training set: ",-log10(pvalue_training_set)))
print(paste0("concordance of training set: ",cindex_training_set[1]))

#
#bring in targeted counts for the test data
#
testCountDF <- readRDS(file.path(normalizePath("."),"derived_phenotype_matrix.rds"))
sel_features <- setdiff(colnames(testCountDF),"name")
testCountDFWM <- inner_join(testCountDF,metaData,by=c("name"))

#
#subset to stim samples
#
stimDF <- testCountDFWM[which(testCountDFWM$stimStatus==1),]
stimDF <- stimDF[,-which(colnames(stimDF)=="stimStatus")]
colnames(stimDF)[which(colnames(stimDF) %in% sel_features)] <- paste0("S_",colnames(stimDF)[which(colnames(stimDF) %in% sel_features)])

#
#subset to unstim samples
#
unstimDF <- testCountDFWM[which(testCountDFWM$stimStatus==0),]
unstimDF <- unstimDF[,-which(colnames(unstimDF)=="stimStatus")]
colnames(unstimDF)[which(colnames(unstimDF) %in% sel_features)] <- paste0("U_",colnames(unstimDF)[which(colnames(unstimDF) %in% sel_features)])

#
#derive features for prediction
#
testDF <- inner_join(stimDF,unstimDF,by=c("subjectID","Status","Survival.Time","inTraining"))
for (cellPop in sel_features) {
    testDF[,paste0("S_",cellPop)] <- (testDF[,paste0("S_",cellPop)]-testDF[,paste0("U_",cellPop)])
}
testDF$Survival.Time <- testDF$Survival.Time-min(min(testDF$Survival.Time),-10)+1

#
#make predictions on test set
#
to_predict_test_set <- testDF[,final_selected_features,drop=FALSE]
prediction_test_set <- predict(rf_training_fit,newdata=to_predict_test_set, na.action="na.impute")

#
#score predictions on test set
#
test_survivalTime <- testDF[,"Survival.Time"]
test_status <- testDF[,"Status"]


fitted_value_test_set <- fitModel(test_survivalTime,test_status,prediction_test_set$predicted)
pvalue_test_set <- summary(fitted_value_test_set)$logtest[3]
cindex_test_set <- summary(fitted_value_test_set)$concordance
print("Analysis results.")
print(paste0("-log10(pvalue) on test set: ",-log10(pvalue_test_set)))
print(paste0("concordance on test set: ",cindex_test_set[1]))



