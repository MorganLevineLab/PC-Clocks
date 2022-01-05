# SPECIFY DIRECTORIES
workingDir <- "[your-path-here/]" #This is where you would like any save data from this script to reside and must
#end in '/'
clocksDir <- "[your-path-here]/PC-Clocks-distribution/"

#Note: You may need to alter the code as needed based on file locations.
#Load packages
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest("dplyr")
library(dplyr)
pkgTest("tibble")
library(tibble)
pkgTest("tidyr")
library(tidyr)


#Designate DNAm and phenotype data
#In datMeth, row names are samples and column names are CpGs.
#In datPheno, rows are samples and columns are phenotypic variables. 
#One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
#Also ensure that the order of datMeth sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled

#Load example data
load(file = paste(clocksDir,"Example_PCClock_Data.RData",sep=""))

#Example code for re-coding females as 1 and males as 0
datPheno$Female <- 2 - as.numeric(as.factor(datPheno$gender))
datPheno$Age <- as.numeric(datPheno$age)
colnames(datPheno)[1] <- "SampleID"
#Check that your samples are in the same order in both methylation and phenotype data
sum(rownames(datMeth) == datPheno$SampleID)


###################
###Optional: calculate original CpG-based clocks (for comparison to PC clocks)
###################
#calculate DNAmAge using online calculator to get GrimAge
datMethOnline <- as.data.frame(t(datMeth)) %>% rownames_to_column("ProbeID")
#Double check that order of datMeth sample IDs match phenotype data samples IDs. If it doesn't match, check that class of datPheno$SampleID is "character"
colnames(datMethOnline)[-1] == datPheno$SampleID
#Filter CpGs to just the relevant ones so that you can upload it to the website (otherwise methylation file is too big)
cgHorvathNew <- read.csv(paste(clocksDir,"cgHorvathNew.csv",sep=""))
subsetCG <- function(dat,cgSet){
  match1 = match(as.matrix(cgSet[,1]),as.matrix(dat[,1]))
  datReduced = dat[match1,]
  datReduced[is.na(match1),1]= as.character(cgSet[is.na(match1),1])
  datReduced
}
datMeth_Horvath_CpGs_New <- subsetCG(datMethOnline,cgHorvathNew)
sum(colnames(datMeth_Horvath_CpGs_New)[-1] == datPheno$SampleID)
write.table(datMeth_Horvath_CpGs_New,file = paste(workingDir,"datMeth_Horvath_New.csv",sep=""), row.names=F, sep="," )
write.table(datPheno,file = paste(workingDir,"datPheno.csv",sep=""),row.names=F,sep="," )
rm(datMeth_Horvath_CpGs_New,cgHorvathNew)
#submit to https://dnamage.genetics.ucla.edu/new
#Select "Normalize Data" and "Advanced Analysis" for results that you use moving forward
#To quickly double-check your data is correctly formatted, you could submit first selecting only "Fast Imputation" and not the other two options
#Then you will get results emailed to you quickly. Then go back and select only "Normalize Data" and "Advanced Analysis".

#Load data from online calculator
DNAmAge <-read.csv(file = paste(workingDir,"datMeth_Horvath_New.output.csv",sep=""))
calcColumns <- c("SampleID","Age","Female",
                "DNAmAge","DNAmAgeSkinBloodClock","DNAmAgeHannum","DNAmPhenoAge","DNAmTL","DNAmGrimAge",
                "DNAmADM","DNAmB2M","DNAmCystatinC","DNAmGDF15","DNAmLeptin","DNAmPAI1","DNAmTIMP1","DNAmPACKYRS")
DNAmAge <- DNAmAge[,calcColumns]
newColnames <- c("SampleID","Age","Female",
                "Horvath1","Horvath2","Hannum","PhenoAge","DNAmTL","GrimAge",
                "DNAmADM","DNAmB2M","DNAmCystatinC","DNAmGDF15","DNAmLeptin","DNAmPAI1","DNAmTIMP1","DNAmPACKYRS")
colnames(DNAmAge) <- newColnames
#if the Horvath calculator added an "X" to the beginning of your SampleIDs:
DNAmAge$SampleID <- sub("X","",DNAmAge$SampleID)

###################
###Calculate PC clocks
###################

load(file = paste(clocksDir,"CalcAllPCClocks.RData",sep=""))

#If needed: Fill in missing CpGs needed for calculation of PCs; use mean values from GSE40279 (Hannum 2013; blood)
datMeth <- as.data.frame(datMeth)
missingCpGs <- CpGs[!(CpGs %in% colnames(datMeth))]
datMeth[,missingCpGs] <- NA
for(i in 1:length(missingCpGs)){
  datMeth[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
}

#Prepare methylation data for calculation of PC Clocks (subset to 78,464 CpGs and perform imputation if needed)
datMeth <- datMeth[,CpGs]
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMeth <- apply(datMeth,2,meanimpute)
#Note: you may substitute another imputation method of your choice (e.g. KNN), but we have not found the method makes a significant difference.

#If you did not calculate the original CpG-based clocks, initialize a data frame for PC clocks
#DNAmAge <- datPheno

#Double check that order of samples is still correct
all(DNAmAge$SampleID == datPheno$SampleID)
  #if FALSE, run sum(DNAmAge$SampleID == datPheno$SampleID) to get how many ARE in order
all(DNAmAge$SampleID == rownames(datMeth))
  #if FALSE, run sum(DNAmAge$SampleID == rownames(datMeth)) to get how many ARE in order

#Calculate PC Clocks
DNAmAge$PCHorvath1 <- anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath1$center) %*% CalcPCHorvath1$rotation %*% CalcPCHorvath1$model + CalcPCHorvath1$intercept)
DNAmAge$PCHorvath2 <- anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath2$center) %*% CalcPCHorvath2$rotation %*% CalcPCHorvath2$model + CalcPCHorvath2$intercept)
DNAmAge$PCHannum <- sweep(as.matrix(datMeth),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept
DNAmAge$PCPhenoAge <- sweep(as.matrix(datMeth),2,CalcPCPhenoAge$center) %*% CalcPCPhenoAge$rotation %*% CalcPCPhenoAge$model + CalcPCPhenoAge$intercept
DNAmAge$PCDNAmTL <- sweep(as.matrix(datMeth),2,CalcPCDNAmTL$center) %*% CalcPCDNAmTL$rotation %*% CalcPCDNAmTL$model + CalcPCDNAmTL$intercept
temp <- cbind(sweep(as.matrix(datMeth),2,CalcPCGrimAge$center) %*% CalcPCGrimAge$rotation,Female = DNAmAge$Female,Age = DNAmAge$Age)
DNAmAge$PCPACKYRS <- temp[,names(CalcPCGrimAge$PCPACKYRS.model)] %*% CalcPCGrimAge$PCPACKYRS.model + CalcPCGrimAge$PCPACKYRS.intercept
DNAmAge$PCADM <- temp[,names(CalcPCGrimAge$PCADM.model)] %*% CalcPCGrimAge$PCADM.model + CalcPCGrimAge$PCADM.intercept
DNAmAge$PCB2M <- temp[,names(CalcPCGrimAge$PCB2M.model)] %*% CalcPCGrimAge$PCB2M.model + CalcPCGrimAge$PCB2M.intercept
DNAmAge$PCCystatinC <- temp[,names(CalcPCGrimAge$PCCystatinC.model)] %*% CalcPCGrimAge$PCCystatinC.model + CalcPCGrimAge$PCCystatinC.intercept
DNAmAge$PCGDF15 <- temp[,names(CalcPCGrimAge$PCGDF15.model)] %*% CalcPCGrimAge$PCGDF15.model + CalcPCGrimAge$PCGDF15.intercept
DNAmAge$PCLeptin <- temp[,names(CalcPCGrimAge$PCLeptin.model)] %*% CalcPCGrimAge$PCLeptin.model + CalcPCGrimAge$PCLeptin.intercept
DNAmAge$PCPAI1 <- temp[,names(CalcPCGrimAge$PCPAI1.model)] %*% CalcPCGrimAge$PCPAI1.model + CalcPCGrimAge$PCPAI1.intercept
DNAmAge$PCTIMP1 <- temp[,names(CalcPCGrimAge$PCTIMP1.model)] %*% CalcPCGrimAge$PCTIMP1.model + CalcPCGrimAge$PCTIMP1.intercept
DNAmAge$PCGrimAge <- as.matrix(DNAmAge[,CalcPCGrimAge$components]) %*% CalcPCGrimAge$PCGrimAge.model + CalcPCGrimAge$PCGrimAge.intercept
rm(CalcPCHorvath1,CalcPCHorvath2,CalcPCHannum,CalcPCPhenoAge,CalcPCDNAmTL,CalcPCGrimAge,temp,imputeMissingCpGs)

#Calculate age residuals if desired. Change clockColumns to the column indices for those clocks you wish to calculate age residuals for
clockColumns = c(4:31)
for (i in clockColumns){
  DNAmAge[,paste0(colnames(DNAmAge)[i],"Resid")] = resid(lm(DNAmAge[,i] ~ DNAmAge$Age))
}
