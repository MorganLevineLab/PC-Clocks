calcPCClocks_Accel <- function(DNAmAge){ 
  #Calculate age residuals
  #DNAmAge is the dataframe where PCClocks were stored
  var = readline(prompt = "Do you want all of the GrimAge Components? (yes/no)")
  if(var == "yes" | var == "Yes"){
    clockColumns = c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCDNAmTL",
                     "PCDNAmTL", "PCPACKYRS", "PCADM", "PCB2M", "PCCystatinC",
                     "PCGDF15", "PCLeptin", "PCPAI1", "PCTIMP1", "PCGrimAge")
  } else {
    clockColumns = c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCDNAmTL", "PCGrimAge")
  }
  #clockColumns = c(9:23)
  for (i in clockColumns){
    DNAmAge[,paste0(i,"Resid")] = resid(lm(DNAmAge[,i] ~ DNAmAge$Age))
    #DNAmAge[,paste0(colnames(DNAmAge)[i],"Resid",sep="")] = resid(lm(as.vector(DNAmAge[,i][[1]]) ~ as.vector(DNAmAge$Age)))
  }
  return(DNAmAge)
}
