# Load the Functions you need to calculate the PCClocks (you need to change the path to the directory 
#       where you installed the code)

clocksDir <- "~/Downloads/PC-Clocks-Beta-Jan2022/PC-Clocks-update-beta/" #where the clocks directory was downloaded to
                                                                         #must end with '/'

source(paste(clocksDir, "run_calcPCClocks.R", sep = ""))
source(paste(clocksDir, "run_calcPCClocks_Accel.R", sep = ""))
load(paste(clocksDir,"CalcAllPCClocks.RData",sep=""))

# Load the file with your pheno data and methylation data in it (Here we used the example data included)
load(paste(clocksDir,"Example_PCClock_Data_final.RData",sep=""))

# IMPORTANT FORMATTING NOTE: If you are not using the example methylation and Pheno data, you will need to have specific
#     formatting. Please ensure that your Methylation dataframe/ matrix is of the methylation beta values and row names
#     are sample names, and column names are CpGs.
#     For the pheno data, ensure that the data frame/ matrix has rows as samples, and columns as whatever phenotype
#     variables you have/ wish. This can also include the original CpG clocks if you used the online Horvath calculator
#     as well. HOWEVER, the pheno matrix MUST have a column named "Age", and a column named "Female" (capital required),
#     especially if you want to calculate GrimAge and its components. Female = 1 is F and Female = 0 is M.
#
#     If you don't have this information, you can also just set it so Females is all 1 (all samples labeled female) and
#     all the same Age. Just know that this won't be accurate for PCGrimAge or components, and that you can't run
#     the acceleration calculations with calcPCClock_Accel.
#
#     The code below is going to ask if you would like to check the order of your datPheno and datMeth samples to ensure
#     they line up. For this to work, you will need to type the column name of datPheno with the names of the samples or 
#     'skip'.

# Get the PC Clocks values and the PC Clock Acceleration values
PCClock_DNAmAge <- calcPCClocks(path_to_PCClocks_directory = clocksDir,
                                     datMeth = datMeth,
                                     datPheno = datPheno)
#in order to calculate Acceleration below, you will need a column called "Age", just as was needed for PCGrimAge.
PCClock_DNAmAge <- calcPCClocks_Accel(PCClock_DNAmAge)


