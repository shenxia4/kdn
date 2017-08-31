# This function is used to run the kernel based neural network on a test data set
# using the weights estimated from the training data set
getTestPredErr <- function(optLambdajVec, optLambdaliMat, optPhi,      # The optimum values for the parameters
                           testBaseKernelList, testInnerKernelName,    # Set the base kernels and inner kernels to be used in the test data
                           testTrait, nSamp = 1e3)
{
  y <- testTrait;
  
  optValList <- CalcDerivLoss(lambdajVec = optLambdajVec, lambdaliMat = optLambdaliMat, phi = optPhi,
                                baseKernelList = testBaseKernelList, innerKernelName = testInnerKernelName,
                                trait = testTrait, nSamp = nSamp);

  
  # Calculate prediction error for test data
  optAMat   <- optValList$Amat;
  optPMat   <- optValList$Pmat
  predErr   <- t(y) %*% optAMat %*% optAMat %*% y;
  y.pred    <- optPMat %*% y;
  predErrL1 <- sum(abs(y - y.pred));
  corY      <- cor(y, y.pred)
  
  optimList <- list(TestPredErr = predErr, TestPredErrL1 = predErrL1, TestCorY = corY);
  return(optimList)
  
}


