#Simulation Setup Codes
simSetUp <- function()
{
	#Load All Necessary sources
  source("./Source/HelperFunctions.R")
  source("./Source/KernelPool.R")
	source("./Source/simFunctionNew.R")
  source("./Source/SingleLayerMultipleKernelNew.R")
  source("./Source/TraitSimulation.R")
  source("./Source/TestDataPredErr.R")



	#Load All Necessary libraries
	library(Matrix)
  library(MASS)


	###################################################################################################
	#Read Data
	###################################################################################################
	ped <- read.table("./data/ped.ped", sep= "\t");
	info <- read.table("./data/info.info", head=T);

	return(list(ped=ped, info=info));
}
