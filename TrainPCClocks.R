library(glmnet)

#To start, you must designate the following:
#CpGs: a vector of CpGs you want to use for your PC clock. Ideally these CpGs should be present in all 
#     of your training and test data sets.
#You can limit the number of CpGs to make calculations faster. 
#For reference, we have found that at least 8,000 CpGs are needed for a reliable predictor of age, and 
#     50,000 CpGs are needed for a reliable predictor of mortality. This number likely varies depending 
#     on your phenotype and on the number of samples.  
#datMethTrain and datMethTest: Training and test DNA methylation data sets, with samples as rows and CpGs as columns
#datPhenoTrain and datPhenoTest: Training and test phenotypic data sets containing the phenotype to be predicted, 
#     with samples in the same order as methylation data.

#Subset methylation data to your set of CpGs
CpGs <- #your-vector-here
datMethTrain <- datMethTrain[,CpGs]
#Impute missing values if needed. You can use a different imputation method of your choice but we have not found 
#     this makes a significant difference.
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain <- apply(datMethTrain,2,meanimpute)

#check order of data is correct
if(all(colnames(datMethTrain) == CpGs)){
  message("CpGs are all in order")
} else(message(paste("Only",sum(colnames(datMethTrain) == CpGs),"CpGs are in order")))

if(all(rownames(datMethTrain) == rownames(datPhenoTrain))){ #may need to change rownames(datPhenoTrain) to a column name
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(datMethTrain) == rownames(datPhenoTrain)),"Samples are in order")))

#Perform PCA and projections. Remove last PC.
PCA = prcomp(datMethTrain,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TestPCData = predict(PCA,datMethTest)[,1:(dim(PCA$x)[2]-1)]

#Select phenotype to be predicted. For example, here we are predicting age.
TrainAge = datPhenoTrain$Age
TestAge = datPhenoTest$Age

#Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainAge, nfolds=10,alpha=0.5, family="gaussian")
fit = glmnet(TrainPCData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)

#Examine full model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.min))

#Examine sparse model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.1se))

#Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCAge <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCAge$model = temp[temp!=0,][-1]
CalcPCAge$intercept = temp[1,1]
CalcPCAge$center = PCA$center
CalcPCAge$rotation = PCA$rotation[,names(CalcPCAge$model)]

#And save your model:
save(CalcPCAge,CpGs,file = "CalcPCAge.RData")

#Then you can calculate your new clock in test data using compressed code, which will now take less time to load and calculate PCs:
load(file = "CalcPCAge.RData")
datMethTest = datMethTest[,CpGs]
PCAge <- sweep(as.matrix(datMethTest),2,CalcPCAge$center) %*% CalcPCAge$rotation %*% CalcPCAge$model + CalcPCAge$intercept

