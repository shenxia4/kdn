# source("HelperFunctions.R")
# source("KernelPool.R")

# Function that calculate the combination of kernel matrices
# HUid is the index for hidden units
CalcBaseKernelComb <- function(lambdaliMat, baseKernelList, HUid)
{
  n <- dim(baseKernelList[[1]])[1];  # number of individuals
  L <- dim(lambdaliMat)[1];              # number of base kernel matrices
  
  IdMat <- diag(rep(1,n));
  
  ResultMatList <- mapply(FUN = '*', baseKernelList, exp(lambdaliMat[, HUid]), SIMPLIFY = FALSE);
  ResultMat <- Reduce(f = '+', ResultMatList) + IdMat;
  
  return(ResultMat)
}


# Function that calculate the combination of kernel matrices in the inner layer
CalcInnerKernelComb <- function(lambdajVec, innerKernelList)
{
  n <- dim(innerKernelList[[1]])[1];    # number of individuals;
  J <- length(lambdajVec);            # number of inner layer kernel matrices
  
  IdMat <- diag(rep(1,n));
  
  ResultMatList <- mapply(FUN = '*', innerKernelList, exp(lambdajVec), SIMPLIFY = FALSE);
  ResultMat <- Reduce(f = '+', ResultMatList) + IdMat
  
  return(ResultMat)
  
}


# Function used to calculate necessary derivatives for gradient descent algorithm
# lambdajVec is the initial vector for the variance components for final prediction;
# lambdaliMat is an L by m matrix containing the initial values for the varaince component for the single layer;
# L is the number of base kernel matrices;
# J is the number of inner layer kernel matrices;
# m is the number of hidden units;
# trait is the phenotype vector;
# baseKernelList is a matrix list containing the L base kernel matrices;
# innerKernelName list the name of all the kernels to be used in the inner layer;
# baseKernelList can have more than L elements and innerKernelName can have more than J elements. If so, the first L and J elements will be used;
CalcDerivMultiKernel <- function(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3)
{
  L <- dim(lambdaliMat)[1]; # number of base kernel matrices
  J <- length(lambdajVec); # number of inner layer matrices
  m <- dim(lambdaliMat)[2]; # number of hidden units;
  n <- length(trait);   # number of individuals
  
  
  # Initialization variables
  IdMat <- diag(rep(1,n));
  Amat <- DevAPhi <- matrix(0, nrow = n, ncol = n);   # used to store results for matrix A
  Pmat <- matrix(0, nrow = n, ncol = n);              # used to store results for matrix P to calculate the predictied y
  DevALambdajList <- InitializeList(m = J, n = n);    # used to store results for derivative of A and A^2 wrt lambda's
  DevALambdaliList <- InitializeList(m = L * m, n = n);
  baseKernelCombList <- list();
  baseKernelCombInvList <- list();
  innerKernelCombList <- list();
  innerKernelCombInvList <- list();
  sampleUMat <- matrix(0, nrow = n, ncol = m);
  
  # Choose the first L and first J kernels in baseKernelList and innerKernelList respectively
  baseKernelList <- baseKernelList[1:L];
  innerKernelList <- list();
  innerKernelName <- innerKernelName[1:J];
  
  # Calculate the covariance matrices for the m hidden units (DOES NOT DEPEND ON U MATRIX!)
  for (i in 1:m)
  {
    baseKernelCombList[[i]] <- CalcBaseKernelComb(lambdaliMat = lambdaliMat, baseKernelList = baseKernelList, HUid = i);
    baseKernelCombInvList[[i]] <- FastInverseMatrix(baseKernelCombList[[i]]);
  }
  
  # Do sampling here to calculate expectations
  for(s in 1:nSamp)
  {
    trMatPhi <- 0;
    
    # Sample U from the multivariate normal distribution
    for(i in 1:m)
    {
      temp <- try(mvrnorm(n = 1, mu = rep(0,n), Sigma = exp(phi) * baseKernelCombList[[i]]));
      if(inherits(temp, 'try-error'))
      {
        print(temp);
      }
      sampleUMat[,i] <- temp;
      trMatPhi <- trMatPhi + baseKernelCombInvList[[i]] %*% tcrossprod(sampleUMat[,i]);
    }
    
    # Obtain the inner layer kernel matrices
    for (j in 1:J)
    {
      innerKernelList[[j]] <- findKernel(innerKernelName[j], geno = sampleUMat);
    }
    
    # Calculate the linear combination of inner layer kernel matrices
    innerKernelCombList[[s]] <- CalcInnerKernelComb(lambdajVec = lambdajVec, innerKernelList = innerKernelList);
    innerKernelCombInvList[[s]] <- FastInverseMatrix(innerKernelCombList[[s]]);
    
    # For matrix A
    Amat <- Amat + (1/nSamp) * innerKernelCombInvList[[s]];
    Pmat <- Pmat + (1/nSamp) * (innerKernelCombList[[s]] - IdMat) %*% innerKernelCombInvList[[s]];
    
    # For derivatives of A with respect to lambda_j
    for (j in 1:J)
    {
      DevALambdajList[[j]] <- DevALambdajList[[j]] - (1/nSamp) * exp(lambdajVec[j])  * innerKernelCombInvList[[s]] %*% 
                                             innerKernelList[[j]] %*%  innerKernelCombInvList[[s]];
    }
    
    # For derivatives of A with respect to lambda_{li}
    for(k in 1:(L*m))
    {
      Lindx <- getIndex(k, L, m)[1];
      mindx <- getIndex(k, L, m)[2];
      
      trMat <- baseKernelCombInvList[[mindx]] %*% baseKernelList[[Lindx]] %*% baseKernelCombInvList[[mindx]] %*% 
                (baseKernelCombList[[mindx]] - exp(-phi) * crossprod(t(sampleUMat[,mindx])));
      trVal <- tr(trMat);
      DevALambdaliList[[k]] <- DevALambdaliList[[k]] - (1/nSamp) * (exp(lambdaliMat[Lindx, mindx]) / 2) * trVal * innerKernelCombInvList[[s]];
    }
    
    
    # For derivative of phi
    DevAPhi <- DevAPhi - (1/(2 * nSamp)) * (m * n - exp(-phi) * tr(trMatPhi)) * innerKernelCombInvList[[s]];
    
  }
  
  returnList <- list(L = L, m = m, n = n, J = J,
                     baseKernelCombList = baseKernelCombList, baseKernelCombInvList = baseKernelCombInvList,
                     innerKernelList = innerKernelList, innerKernelCombList = innerKernelCombList, innerKernelCombInvList = innerKernelCombInvList,
                     Amat = Amat, DevALambdajList = DevALambdajList, DevALambdaliList = DevALambdaliList, DevAPhi = DevAPhi,
                     Pmat = Pmat);
  
  return(returnList)
  
}


# Function for calculating derivatives with respect to the loss function
CalcDerivLoss <- function(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3)
{
  matDerivList <- CalcDerivMultiKernel(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp);
  L <- matDerivList$L; J <- matDerivList$J;
  m <- matDerivList$m; n <- matDerivList$n;
  y <- trait;
  
  Amat <- matDerivList$Amat;
  DevALambdajList <- matDerivList$DevALambdajList;
  DevALambdaliList <- matDerivList$DevALambdaliList;
  DevAPhi <- matDerivList$DevAPhi;
  
  derivLoss <- rep(0, 1 + J + L * m);
  
  for(i in 1:length(derivLoss))
  {
    if(i <= J)
    {
      tempMat <- DevALambdajList[[i]] %*% Amat + Amat %*% DevALambdajList[[i]];
      derivLoss[i] <- t(y) %*% tempMat %*% y;
    }else
      if(i <= J + L*m)
      {
        k <- i - J;
        tempMat <- DevALambdaliList[[k]] %*% Amat + Amat %*% DevALambdaliList[[k]];
        derivLoss[i] <- t(y) %*% tempMat %*% y;
      }else
      {
        tempMat <- DevAPhi %*% Amat + Amat %*% DevAPhi;
        derivLoss[i] <- t(y) %*% tempMat %*% y;
      }
  }
  
  lossDerivList <- list(derivLoss = derivLoss);
  returnList <- c(matDerivList, lossDerivList);
  return(returnList)
}


# Function for gradient descent
# niter is the maximum iterations for gradient descent
# lr is a fixed learning rate used in the option specified and the first iteration in the option of B-Bmethod
GradDesc <- function(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3, niter = 100, lr = 0.001, tol = 1e-5)
{
  derivList <- CalcDerivLoss(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp);
  L <- derivList$L; J <- derivList$J;
  m <- derivList$m; n <- derivList$n;
  y <- trait;
  
  # Gradients
  lossDeriv <- derivList$derivLoss;
  
  # Combine all parameters in a vector
  paraVec <- c(c(lambdajVec), c(t(lambdaliMat)), phi);
  
  # Matrix used to store the iteration result
  paraMat <- matrix(0, nrow = niter, ncol = length(paraVec));
  paraMat[1,] <- paraVec;
  
  for(i in 2:niter)
  {
    cat("Iteration", i, "\n");
    if(i == 2)
    {
      lr <- getLearningRate(lr = lr, type = "Specified");
      paraMat[i,] <- paraMat[i-1,] - lr * lossDeriv;
      devFxn <- lossDeriv;
    }else
    {
      xn <- paraMat[i-1,]; xn1 <- paraMat[i-2,];
      devFxn1 <- devFxn;
      
      lambdajVecIter <- paraMat[i-1, 1:J];
      lambdaliMatIter <- matrix(paraMat[i-1, (J+1):(J+L*m)], nrow = L, byrow = T);
      phiIter <- paraMat[i-1, length(paraVec)];
      
      derivListIter <- CalcDerivLoss(lambdajVec = lambdajVecIter, lambdaliMat = lambdaliMatIter, phi = phiIter,
                                     baseKernelList = baseKernelList, innerKernelName = innerKernelName, trait = trait, nSamp = nSamp);
      lossDerivInter <- derivListIter$derivLoss;
      devFxn <- lossDerivInter;
      
      lr <- getLearningRate(xn, xn1, devFxn, devFxn1, type = "B-BMethod");
      
      paraMat[i,] <- paraMat[i-1,] - lr * lossDerivInter;
    }
    
    cat("Learning Rate = ", lr, "\n");
    cat("lambdaj = [", paraMat[i, 1:J], "], lambdali = [", paraMat[i, (J+1):(J+L*m)], "], phi = ", paraMat[i, length(paraVec)], "\n");
    
    if(sum(abs(paraMat[i,] - paraMat[i-1,])) < tol)
    {break;}
  }
  
  gradDescList <- list(lambdajVec = paraMat[i, 1:J], lambdaliMat = matrix(paraMat[i, (J+1):(J+L*m)], nrow = L, byrow = T), phi = paraMat[i, length(paraVec)]);
  returnList <- c(derivList, gradDescList);
  return(returnList)
}


# Function for prediction error
CalcPredErr <- function(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3, niter = 100, lr = 0.001, tol = 1e-5)
{
  y <- trait;
  gradDescList <- GradDesc(lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp, niter, lr, tol);
  
  lambdajVecOpt <- gradDescList$lambdajVec;
  lambdaliMatOpt <- gradDescList$lambdaliMat;
  phiOpt <- gradDescList$phi;
  
  optValList <- CalcDerivLoss(lambdajVec = lambdajVecOpt, lambdaliMat = lambdaliMatOpt, phi = phiOpt,
                              baseKernelList = baseKernelList, innerKernelName = innerKernelName, trait = trait, nSamp = nSamp);
  
  # Calculate prediction error
  optAMat   <- optValList$Amat;
  optPMat   <- optValList$Pmat
  predErr   <- t(y) %*% optAMat %*% optAMat %*% y;
  y.pred    <- optPMat %*% y;
  predErrL1 <- sum(abs(y - y.pred));
  corY      <- cor(y, y.pred)
  
  optimList <- list(lambdajVec = lambdajVecOpt, lambdaliMat = lambdaliMatOpt, phi = phiOpt, 
                    predErr = predErr, predErrL1 = predErrL1, corY = corY);
  returnList <- c(optValList, optimList);
  return(returnList)
}
