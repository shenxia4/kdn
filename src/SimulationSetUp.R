## Simulation Setup Codes
simSetUp <- function()
{
    ## Load All Necessary sources
    source("src/HelperFunctions.R")
    source("src/KernelPool.R")
    source("src/simFunctionNew.R")
    source("src/SingleLayerMultipleKernelNew.R")
    source("src/TraitSimulation.R")

    ## Load All Necessary libraries
    library(Matrix)
    library(MASS)


    ## ------------  Read Data  ------------ ##
    ped <- read.table("dat/ped.ped", sep= "\t")
    nfo <- read.table("dat/nfo.nfo", head=T)

    return(list(ped=ped, nfo=nfo))
}
