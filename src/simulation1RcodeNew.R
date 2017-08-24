## Simulation Setup 1 Codes
source("src/HelperFunctions.R")
source("src/KernelPool.R")
source("src/simFunctionNew.R")
source("src/sgd.R")
source("src/TraitSimulation.R")

#' the main evaluation entrance
#' genotype data is loaded from 'dat/gno.rds'
#' 
#' @param N   integer, samples to be generated
#' @param ffr numeric, the proportion of functional variants, aka Function
#' FRequency
#' @param lgs integer, number of genomic variants to chop from raw data.
#'
#' @param itr integer, number of simulation to run, aka the outer most iteration.
#' @param nep integer, number of allowed epoches, aks the intermediate iteration.
#' @param mvn integer, number of multi variate normal to be sampled for each inner
#' unit, aka the inner most iteration.
main <- function(N=200, ffr=.5, itr=1, ngv=3e4, nep=100, mvn=20)
{
    dat <- readRDS('dat/gno.rds')
    ped <- dat$ped
    nfo <- dat$nfo

    tag <- 2

    ## trait simulated based on identity matrix, the covariance matrix of random effect
    traitSigma <- diag(rep(1,N))
    traitFun <- "identity"
    sigmaR <- 1
    phi <- 1
    order <- 1

    ## Parameters for gradient descent algorithm
    fromBaseKernel <- c("CAR", "identity")
    ## innerKernel <- c("identity", "product")
    innerKernel <- c("product", "product", "identity", "product")
    lambdajVec <- c(1,1)
    lambdaliMat <- matrix(c(0,1,2,0,1,2), nrow = 2, byrow = T)
    tol <- 1e-10

    Rprof()
    time1 <- proc.time()
    pred1 <- getPred(
        tag = tag, func.frq = ffr, nS = itr, nSeq = ngv, nN = N, info = nfo, ped = ped,
        traitMu = rep(0, N), traitSigma = traitSigma, traitFun = traitFun,
        sigmaR = sigmaR, phi = phi, order = order,
        fromBaseKernel = fromBaseKernel, UserDef = NULL, 
        lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
        innerKernel = innerKernel,
        nSamp = mvn, niter = nep, tol = tol)
    proc.time()-time1
    Rprof(NULL)

    pred0 <- readRDS('dat/s525c.rds')
    all.equal(pred0, pred1)

    list(ref=pred0, out=pred1)
}
