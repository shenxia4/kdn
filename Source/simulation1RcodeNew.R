## Simulation Setup 1 Codes
## setwd("C:/Users/xshen/Desktop/Reasearch Projects/Deep Learning Project/Simulation/Simulation1/DeepLearningSimulation1")
source("./Source/SimulationSetUp.R")

dataList <- simSetUp();
ped <- dataList$ped;
info <- dataList$info;


tag <- 1;
func.frq <- 0.5;         #number of functional alleles

nS <- 100;     #number of simulations
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

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,1,2,3,4,0,1,2,3,4), nrow = 2, byrow = T);
nSamp <- 200;
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

#####################################################################################################################################################################
## Simulation 2: On data generated using CAR kernel
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "CAR";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,1,2,0,1,2), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;

time1 <- proc.time()
pred2 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1

##########################################################################################################################################################
## Simulation 3: On data generated using product kernel
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "product";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,1,2,0,1,2), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;

time1 <- proc.time()
pred3 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1


##########################################################################################################################################################
## Simulation 4: On data generated using product kernel with base kernel having product kernel
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "product";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity", "product");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(rep(0,9), nrow = 3, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;

time1 <- proc.time()
pred4 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1


##########################################################################################################################################################
## Simulation 5: Nonlinear relationship
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "CAR";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "polynomial";
sigmaR <- 1;
phi <- 1;
order <- 2;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,0,0,0,0,0), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;


time1 <-proc.time()
pred5 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1

##########################################################################################################################################################
## Simulation 6: Nonlinear relationship
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "CAR";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "sin";
sigmaR <- 1;
phi <- 1;
order <- 2;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(c(0,0,0,0,0,0), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;


time1 <-proc.time()
pred6 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1



##########################################################################################################################################################
## Simulation 7: On data generated using Gaussian kernel without base kernel having Gaussian kernel
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "Gaussian";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("CAR", "identity");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(rep(0,6), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;

time1 <- proc.time()
pred7 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1


##########################################################################################################################################################
## Simulation 8: On data generated using Gaussian kernel with base kernel having Gaussian kernel
# Parameters for generating traits
traitMu <- rep(0, nN); 
traitSigma <- "Gaussian";  # trait simulated based on identity matrix, the covariance matrix of random effect
traitFun <- "identity";
sigmaR <- 1;
phi <- 1;
order <- 1;

# Parameters for gradient descent algorithm
fromBaseKernel <- c("identity", "Gaussian");
innerKernelName <- c("identity", "product");
lambdajVec <- c(1,1);
lambdaliMat <- matrix(rep(0,6), nrow = 2, byrow = T);
nSamp <- 100;
niter <- 100;
tol <- 1e-5;

time1 <- proc.time()
pred8 <- getPred(tag = tag, func.frq = func.frq, nS = nS, nSeq = nSeq, nN = nN, info = info, ped = ped,
                 traitMu = traitMu, traitSigma = traitSigma, traitFun = traitFun,
                 sigmaR = sigmaR, phi = phi, order = order,
                 fromBaseKernel = fromBaseKernel, UserDef = NULL, 
                 lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, varPhi = 1,
                 innerKernelName = innerKernelName,
                 nSamp = nSamp, niter = niter, tol = tol, test = test);
proc.time()-time1


