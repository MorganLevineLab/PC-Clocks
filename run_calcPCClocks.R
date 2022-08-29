calcPCClocks <- function(path_to_PCClocks_directory, datMeth, datPheno){
  #path_to_PCClocks should end with a "/"
  #datMeth is a matrix of methylation Beta values, where row names are samples, and 
  #   column names are CpGs
  #datPheno has rows as samples and columns as phenotype variables. This can also
  #   include the original clocks if you used the Horvath online calculator as well.
  #   It MUST include a column named "Age" and a column named "Female"
  
  if(!("Age" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Age")
  }
  if(!("Female" %in% variable.names(datPheno))){
    stop("Error: datPheno must have a column named Female")
  }
  if(sum(startsWith(colnames(datMeth),"cg")) == 0){
    warning("Warning: It looks like you may need to format datMeth using t(datMeth) to get samples as rows!")
  }


  #Note: this code assumes all your files are in one working directory. Alter the code as needed based on file locations.
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
  

  #In datPheno, rows are samples and columns are phenotypic variables. 
  #One of the phenotypic variables must be "Age", and another one "Female" (coded as Female = 1, Male = 0; should be a numeric variable as this will be included in PCGrimAge calculation)
  #Also ensure that the order of datMeth sample IDs matches your phenotype data sample IDs, otherwise your data will be scrambled
  
  load(file = paste(path_to_PCClocks_directory,"CalcAllPCClocks.RData", sep = ""))
  
  message("PCClocks Data successfully loaded")
  
  #If needed: Fill in missing CpGs needed for calculation of PCs; use mean values from GSE40279 (Hannum 2013; blood)- note that for other tissues you might prefer to use a different one
  datMeth <- as.data.frame(datMeth)
  if(length(c(CpGs[!(CpGs %in% colnames(datMeth))],CpGs[apply(datMeth[,colnames(datMeth) %in% CpGs], 2, function(x)all(is.na(x)))])) == 0){
    message("No CpGs were NA for all samples")
  } else{
    missingCpGs <- c(CpGs[!(CpGs %in% colnames(datMeth))])
    datMeth[,missingCpGs] <- NA
    datMeth = datMeth[,CpGs]
    missingCpGs <- CpGs[apply(datMeth[,CpGs], 2, function(x)all(is.na(x)))]
    for(i in 1:length(missingCpGs)){
      datMeth[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
    }
    message("Any missing CpGs successfully filled in (see function for more details)")
  }
                              
  #Prepare methylation data for calculation of PC Clocks (subset to 78,464 CpGs and perform imputation if needed)
  datMeth <- datMeth[,CpGs]
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  datMeth <- apply(datMeth,2,meanimpute)
  #Note: you may substitute another imputation method of your choice (e.g. KNN), but we have not found the method makes a significant difference.
  message("Mean imputation successfully completed for any missing CpG values")
  
  #Initialize a data frame for PC clocks
  DNAmAge <- datPheno
  
  var = readline(prompt = "To check whether datMeth and datPheno match up, type the column name in datPheno with sample names (or type skip):")
  if(var != "skip"){
    if(sum(DNAmAge[,var] == rownames(datMeth)) != dim(DNAmAge[,var])[1]){
      warning("Warning: It would appear that datPheno and datMeth do not have matching sample order! Check your inputs!")
    } else message("datPheno and datMeth sample order verified to match!")
  }
  
  message("Calculating PC Clocks now")
  
  #Calculate PC Clocks
  DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath1$center) %*% CalcPCHorvath1$rotation %*% CalcPCHorvath1$model + CalcPCHorvath1$intercept))
  DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath2$center) %*% CalcPCHorvath2$rotation %*% CalcPCHorvath2$model + CalcPCHorvath2$intercept))
  DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept)
  DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCPhenoAge$center) %*% CalcPCPhenoAge$rotation %*% CalcPCPhenoAge$model + CalcPCPhenoAge$intercept)
  DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCDNAmTL$center) %*% CalcPCDNAmTL$rotation %*% CalcPCDNAmTL$model + CalcPCDNAmTL$intercept)
  temp <- cbind(sweep(as.matrix(datMeth),2,CalcPCGrimAge$center) %*% CalcPCGrimAge$rotation,Female = DNAmAge$Female,Age = DNAmAge$Age)
  DNAmAge$PCPACKYRS <- as.numeric(temp[,names(CalcPCGrimAge$PCPACKYRS.model)] %*% CalcPCGrimAge$PCPACKYRS.model + CalcPCGrimAge$PCPACKYRS.intercept)
  DNAmAge$PCADM <- as.numeric(temp[,names(CalcPCGrimAge$PCADM.model)] %*% CalcPCGrimAge$PCADM.model + CalcPCGrimAge$PCADM.intercept)
  DNAmAge$PCB2M <- as.numeric(temp[,names(CalcPCGrimAge$PCB2M.model)] %*% CalcPCGrimAge$PCB2M.model + CalcPCGrimAge$PCB2M.intercept)
  DNAmAge$PCCystatinC <- as.numeric(temp[,names(CalcPCGrimAge$PCCystatinC.model)] %*% CalcPCGrimAge$PCCystatinC.model + CalcPCGrimAge$PCCystatinC.intercept)
  DNAmAge$PCGDF15 <- as.numeric(temp[,names(CalcPCGrimAge$PCGDF15.model)] %*% CalcPCGrimAge$PCGDF15.model + CalcPCGrimAge$PCGDF15.intercept)
  DNAmAge$PCLeptin <- as.numeric(temp[,names(CalcPCGrimAge$PCLeptin.model)] %*% CalcPCGrimAge$PCLeptin.model + CalcPCGrimAge$PCLeptin.intercept)
  DNAmAge$PCPAI1 <- as.numeric(temp[,names(CalcPCGrimAge$PCPAI1.model)] %*% CalcPCGrimAge$PCPAI1.model + CalcPCGrimAge$PCPAI1.intercept)
  DNAmAge$PCTIMP1 <- as.numeric(temp[,names(CalcPCGrimAge$PCTIMP1.model)] %*% CalcPCGrimAge$PCTIMP1.model + CalcPCGrimAge$PCTIMP1.intercept)
  DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,CalcPCGrimAge$components]) %*% CalcPCGrimAge$PCGrimAge.model + CalcPCGrimAge$PCGrimAge.intercept)
  rm(CalcPCHorvath1,CalcPCHorvath2,CalcPCHannum,CalcPCPhenoAge,CalcPCDNAmTL,CalcPCGrimAge,temp,imputeMissingCpGs)
  
  message("PC Clocks successfully calculated!")
  
  return(DNAmAge)
}
  
