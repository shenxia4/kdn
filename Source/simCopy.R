## Simulation Setup 1 Codes
source("./Source/SimulationSetUp.R")

dataList <- simSetUp();
ped <- dataList$ped;
info <- dataList$info;


tag <- 1;
func.frq <- 0.5;         #number of functional alleles

nS <- 2;     #number of simulations
nSeq <- 3e4; #length of the genomic region

nN<-100; #or 100, 1000  #number of subjects

test <- 1;


###################################################################################################################################################################
## Simulation 1: On data generated using identity matrix.
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "identity";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

set.seed(100)

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,1,2,3,4,0,1,2,3,4), nrow = 2, byrow = T);
nSamp <- 10;
niter <- 100;
tol <- 1e-4;

time1 <- proc.time()
pred1 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1
