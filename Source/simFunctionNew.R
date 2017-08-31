# Function for getting weights of rare variants
getWeights <- function(tag, geno, maf)
{
  if(tag==1)
  {
    wt<-1;
  }
  if(tag == 2)
  {
    wt<-dbeta(maf,1,25)^2 ;
  }
  if(tag == 3)
  {
    wt<-1 / (maf * (1 - maf));
  } 
  if(tag == 0 || tag == 4)
  {
    wt<--log10(maf);
  }
  
  return(wt);
}

#get best predictor under LMM
#D is the covariance matrix for random effect;
#Z is the design matrix for random effect;
#R is the covariance matrix for error;
getBPLMM <- function(y, SigmaName, Z, geno, maf=NULL, weights = 1, mu.y = 0, mu.rand = 0, R = 0)
{
  D <- findKernel(SigmaName, geno = geno, maf = maf, weights = weights);
  C <- Z %*% D;
  V <- Z %*% D %*% t(Z) + R;
  
  BP <- mu.rand + C %*% solve(V) %*% (y - mu.y);
  
  return(BP)
}

# Function for getting covariance matrix (Kernels)
# fromBaseKernel is a string vector indicating the kernel to choose from the base kenels.
getKernelList <- function(fromBaseKernel , UserDef = NULL, geno, maf = NULL, weights = 1)
{
  kernelList <- list();
  baseKernelName <- c('CAR', 'identity', 'product');
  
  TrueIndx <- baseKernelName %in% fromBaseKernel;
  
  for(i in 1:sum(TrueIndx))
  {
    if(isTRUE(TrueIndx[i]))
    {
      kernelName <- baseKernelName[match(fromBaseKernel[i], baseKernelName)];
      kernelList[[kernelName]] <- findKernel(kernelName, geno, maf, weights);
    }
  }
  
  
  
  ## append user defined kernels into the kernelList
  if(! is.null(UserDef))
  {
    newLength <- length(UserDef);
    for (i in 1:newLength)
    {
      kernelList[[i + length(kernelList)]] <- UserDef[[i]];
    }
  }
  
  return(kernelList)
}


# Main function for simulation
# tag is to indicate the weight used for rare variants;
# func.frq is the frequency of functional variants;
# nS is the number of simulations to be run;
# nN is the total number of individuals;
# test is a tag used to indicate whether to divide the 
# whole data set into training data and test data;
# If test == 1, 
# The first half of nN will be used as training data set;
# The second half of nN will be used as test data set;
getPred <- function(tag, func.frq, nS, nSeq, nN, info, ped, # these variables are used to generate SNP data
                    traitMu, traitSigma, traitFun, 
                    sigmaR, phi, order,                     # these variables are used to generate trait
                    fromBaseKernel, UserDef = NULL,         # these variables are used to generate base kernel matrix list                 
                    lambdajVec, lambdaliMat, varPhi,
                    innerKernelName, 
                    nSamp = 1e3, niter = 100, tol = 1e-5,   # these variables are used to assign initial values to the parameters
                    test = 0)                               # Indicate whether to separate the data into training and testing
{
  predErr <- rep(0, nS);
  predErr.L1 <- rep(0, nS);
  predErrTest <- matrix(0, nS, 2);
  predErrTest.L1 <- matrix(0, nS, 2);
  predErrTest.CorY <- matrix(0, nS, 2);
  predErr.BLUPK <- rep(0, nS);
  predErr.BLUPI <- rep(0, nS);
  predErr.BLUPK.L1 <- rep(0, nS);
  predErr.BLUPI.L1 <- rep(0, nS);
  corY.learn <- rep(0, nS);
  corY.BLUPK <- rep(0, nS);
  corY.BLUPI <- rep(0, nS);
  
  IdMat <- diag(rep(1, nN));
  
  seg.pos <- runif(nS, 0, 1e6-nSeq);
  n.rare <- NULL;
  
  pb <- txtProgressBar(min = 0, max = nS, style = 3);
  cat("\n");
  
  for(i in 1:nS)
  {
    idx <- (info$pos > seg.pos[i] & info$pos < seg.pos[i] + nSeq);	#get a genomic region
    maf <- info$maf[idx];
    
    if(test == 0)
    {
      smp <- sample(1:nrow(ped), nN);						#sample nN individuals
      
      geno <- ped[smp, idx];								#obtain their genetic information within the region
      
      vrt <- apply(geno, 2, variant);
      geno <- as.matrix(geno[,vrt]);
      maf <- maf[vrt];
      n.rare <- c(n.rare, ncol(geno));
      
      
      idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq);
      geno.true<-as.matrix(geno[,idx.true]);
      maf.true<-maf[idx.true];
      
      weights <- getWeights(tag = tag, geno = geno.true, maf = maf.true);
      
      
      # Generate trait list
      traitList <- traitSim(mu = traitMu, Sigma = traitSigma , geno = geno.true, maf = maf.true, weights = weights, 
                            sigmaR = sigmaR, phi = phi, order = order, fun = traitFun);
      trait <- traitList$trait;
      
      # Generate base kernel matrix list
      KList <- getKernelList(fromBaseKernel = fromBaseKernel, UserDef = UserDef, geno = geno.true, maf = maf.true, weights = weights);
      
      # Prediction error calculated here
      predList <- CalcPredErr(lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, phi = varPhi, baseKernelList = KList, 
                              innerKernelName = innerKernelName, trait = trait, nSamp = nSamp, niter = niter, tol = tol);
      predErr[i] <- predList$predErr;
      predErr.L1[i] <- predList$predErrL1;
      corY.learn[i] <- predList$corY;
      
      predVal.BLUPK <- getBPLMM(y = trait, SigmaName = traitSigma, Z = IdMat, geno = geno.true, maf = maf.true, weights = weights, R = phi * IdMat);
      predVal.BLUPI <- getBPLMM(y = trait, SigmaName = "identity", Z = IdMat, geno = geno.true, maf = maf.true, weights = weights, R = phi * IdMat);
      predErr.BLUPK[i] <- t(trait - predVal.BLUPK) %*% (trait - predVal.BLUPK);
      predErr.BLUPI[i] <- t(trait - predVal.BLUPI) %*% (trait - predVal.BLUPI);
      predErr.BLUPK.L1[i] <- sum(abs(trait - predVal.BLUPK));
      predErr.BLUPI.L1[i] <- sum(abs(trait - predVal.BLUPI));
      corY.BLUPK[i] <- cor(trait, predVal.BLUPK);
      corY.BLUPI[i] <- cor(trait, predVal.BLUPI);
      
    }else
    {
      smp       <- sample(1:nrow(ped), (2*nN));						# sample nN individuals
      smp.train <- smp[1:nN];                     # training sample
      smp.test  <- smp[(nN+1):(2*nN)];                # test sample
      
      geno.train <- ped[smp.train, idx];							# obtain their genetic information within the region
      geno.test  <- ped[smp.test, idx];
      
      vrt.train <- apply(geno.train, 2, variant);
      vrt.test  <- apply(geno.test, 2, variant);
      
      geno.train <- as.matrix(geno.train[ , vrt.train]);
      geno.test  <- as.matrix(geno.test[ , vrt.test]);
      
      maf.train <- maf[vrt.train];
      maf.test  <- maf[vrt.test];
      
      n.rare.train <- c(n.rare, ncol(geno.train));
      n.rare.test  <- c(n.rare, ncol(geno.test));
      
      idx.true.train <- sample(1:ncol(geno.train), ncol(geno.train) * func.frq);
      idx.true.test  <- sample(1:ncol(geno.test), ncol(geno.test) * func.frq);
      
      geno.true.train <- as.matrix(geno.train[, idx.true.train]);
      geno.true.test  <- as.matrix(geno.test[, idx.true.test]);
      
      maf.true.train <- maf[idx.true.train];
      maf.true.test  <- maf[idx.true.test];
      
      weights.train <- getWeights(tag = tag, geno = geno.true.train, maf = maf.true.train);
      weights.test  <- getWeights(tag = tag, geno = geno.true.test, maf = maf.true.test);
      
      
      # Generate trait list
      traitList.train <- traitSim(mu = traitMu, Sigma = traitSigma , geno = geno.true.train, maf = maf.true.train, weights = weights.train, 
                            sigmaR = sigmaR, phi = phi, order = order, fun = traitFun);
      traitList.test <- traitSim(mu = traitMu, Sigma = traitSigma , geno = geno.true.test, maf = maf.true.test, weights = weights.test, 
                                  sigmaR = sigmaR, phi = phi, order = order, fun = traitFun);
      trait.train <- traitList.train$trait;
      trait.test  <- traitList.test$trait; 
     
      
      # Generate base kernel matrix list
      KList.train <- getKernelList(fromBaseKernel = fromBaseKernel, UserDef = UserDef, 
                                   geno = geno.true.train, maf = maf.true.train, weights = weights.train);
      KList.test  <- getKernelList(fromBaseKernel = fromBaseKernel, UserDef = UserDef, 
                                   geno = geno.true.test, maf = maf.true.test, weights = weights.test);
      
      # Prediction error calculated here
      # Train the model
      trainList <- CalcPredErr(lambdajVec = lambdajVec, lambdaliMat = lambdaliMat, phi = varPhi, baseKernelList = KList.train, 
                              innerKernelName = innerKernelName, trait = trait.train, nSamp = nSamp, niter = niter, tol = tol);
      optLambdajVec <- trainList$lambdajVec;
      optLambdaliMat <- trainList$lambdaliMat;
      optPhi <- trainList$phi;
      trainErr <- trainList$predErr;
      trainErrL1 <- trainList$predErrL1;
      trainCorY <- trainList$corY;
      
      # Test the model
      testList <- getTestPredErr(optLambdajVec = optLambdajVec, optLambdaliMat = optLambdaliMat, optPhi = optPhi,
                                 testBaseKernelList = KList.test, testInnerKernelName = innerKernelName,
                                 testTrait = trait.test, nSamp = nSamp);
      testErr <- testList$TestPredErr;
      testErrL1 <- testList$TestPredErrL1;
      testCorY <- testList$TestCorY;
      
      predErrTest[i,] <- c(trainErr, testErr);
      predErrTest.L1[i,] <- c(trainErrL1, testErrL1);
      predErrTest.CorY[i,] <- c(trainCorY, testCorY);
      
      predVal.BLUPK <- getBPLMM(y = trait.test, SigmaName = traitSigma, Z = IdMat, geno = geno.true.test, maf = maf.true.test, weights = weights.test, R = phi * IdMat);
      predVal.BLUPI <- getBPLMM(y = trait.test, SigmaName = "identity", Z = IdMat, geno = geno.true.test, maf = maf.true.test, weights = weights.test, R = phi * IdMat);
      predErr.BLUPK[i] <- t(trait.test - predVal.BLUPK) %*% (trait.test - predVal.BLUPK);
      predErr.BLUPI[i] <- t(trait.test - predVal.BLUPI) %*% (trait.test - predVal.BLUPI);
      predErr.BLUPK.L1[i] <- sum(abs(trait.test - predVal.BLUPK));
      predErr.BLUPI.L1[i] <- sum(abs(trait.test - predVal.BLUPI));
      corY.BLUPK[i] <- cor(trait.test, predVal.BLUPK);
      corY.BLUPI[i] <- cor(trait.test, predVal.BLUPI);
      
    }
    
    
    setTxtProgressBar(pb, i);
    cat("\n");
  }
  
  close(pb)
  
  
  if(test == 0)
  {
    result <- matrix(c(mean(predErr), sd(predErr), mean(predErr.L1), sd(predErr.L1), mean(corY.learn), sd(corY.learn),
                       mean(predErr.BLUPK), sd(predErr.BLUPK), mean(predErr.BLUPK.L1), sd(predErr.BLUPK.L1), mean(corY.BLUPK), sd(corY.BLUPK),
                       mean(predErr.BLUPI), sd(predErr.BLUPI), mean(predErr.BLUPI.L1), sd(predErr.BLUPI.L1), mean(corY.BLUPI), sd(corY.BLUPI)), 
                     nrow = 1);
    colnames(result) <- c("Pred.Err.Mean", "Pred.Err.SD", "Pred.Err.L1.Mean", "Pred.Err.L1.SD", "Pred.CorL.Mean", "Pred.CorL.SD",
                          "Pred.Err.Ker.Mean", "Pred.Err.Ker.SD", "Pred.Err.Ker.L1.Mean", "Pred.Err.Ker.L1.SD", "Pred.Err.Ker.Cor.Mean", "Pred.Err.Ker.Cor.SD",
                          "Pred.Err.Id.Mean", "Pred.Err.Id.SD", "Pred.Err.Id.L1.Mean", "Pred.Err.Id.L1.SD", "Pred.Err.Id.Cor.Mean", "Pred.Err.Id.Cor.SD");
  }else
  {
    ErrMean <- apply(predErrTest, 2, mean);
    ErrSD   <- apply(predErrTest, 2, sd);
    L1ErrMean <- apply(predErrTest.L1, 2, mean);
    L1ErrSD   <- apply(predErrTest.L1, 2, sd);
    CorMean <- apply(predErrTest.CorY, 2, mean);
    CorSD   <- apply(predErrTest.CorY, 2, sd);
    result  <- matrix(c(ErrMean[1], ErrSD[1], ErrMean[2], ErrSD[2], 
                        L1ErrMean[1], L1ErrSD[1], L1ErrMean[2], L1ErrSD[2],
                        CorMean[1], CorSD[1], CorMean[2], CorSD[2],
                        mean(predErr.BLUPK), sd(predErr.BLUPK), mean(predErr.BLUPK.L1), sd(predErr.BLUPK.L1), mean(corY.BLUPK), sd(corY.BLUPK),
                        mean(predErr.BLUPI), sd(predErr.BLUPI), mean(predErr.BLUPI.L1), sd(predErr.BLUPI.L1), mean(corY.BLUPI), sd(corY.BLUPI)),
                      nrow = 1);
    colnames(result) <- c("Train.Err.Mean", "Train.Err.SD", "Test.Err.Mean", "Test.Err.SD", 
                          "Train.L1.Err.Mean", "Train.L1.Err.SD", "Test.L1.Err.Mean", "Test.L1.Err.SD",
                          "Train.Cor.Mean", "Train.Cor.SD", "Test.Cor.Mean", "Test.Cor.SD", 
                          "Pred.Err.Ker.Mean", "Pred.Err.Ker.SD", "Pred.Err.Ker.L1.Mean", "Pred.Err.Ker.L1.SD", "Pred.Err.Ker.Cor.Mean", "Pred.Err.Ker.Cor.SD",
                          "Pred.Err.Id.Mean", "Pred.Err.Id.SD", "Pred.Err.Id.L1.Mean", "Pred.Err.Id.L1.SD", "Pred.Err.Id.Cor.Mean", "Pred.Err.Id.Cor.SD");
  }
  
  return(result)
  
}