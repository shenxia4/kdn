#Simulation Setup 1 Codes
source("src/SimulationSetUp.R")
source("src/benchmark.R")

main <- function()
{
    dataList <- simSetUp()
    ped <- dataList$ped
    nfo <- dataList$nfo

    tag <- 2
    func.frq <- 0.5                     # number of functional alleles

    nS <- 2                             # number of simulations
    nSeq <- 3e4                         # length of the genomic region

    nN<-100                             # number of subjects

                                        # Parameters for generating traits
    traitMu <- rep(0, nN)
    ## trait simulated based on identity matrix, the covariance matrix of random effect
    traitSigma <- diag(rep(1,nN))
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
    nSamp <- 10
    niter <- 20
    tol <- 1e-5

    Rprof()
    time1 <- proc.time()
    pred1 <- getPred(
        tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = nfo, ped = ped,
        traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
        sigmaR = sigmaR, phi = phi, order = order,
        fromBaseKernel = fromBaseKernel, UserDef = NULL, 
        lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
        innerKernel = innerKernel,
        nSamp = nSamp, niter = niter, tol = tol)
    proc.time()-time1
    Rprof(NULL)

    pred0 <- readRDS('dat/s525c.rds')
    all.equal(pred0, pred1)

    list(ref=pred0, out=pred1)
}
