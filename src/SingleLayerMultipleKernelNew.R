## source("HelperFunctions.R")
## source("KernelPool.R")

## Function that calculate the combination of kernel matrices
## HUid is the index for hidden units
CalcBaseKernelComb <- function(lambdaliMat, baseKernelList, HUid)
{
    n <- dim(baseKernelList[[1]])[1]; # number of individuals
    L <- dim(lambdaliMat)[1];         # number of base kernel matrices
    
    IdMat <- diag(rep(1,n));
    
    ResultMatList <- mapply(FUN = '*', baseKernelList, exp(lambdaliMat[, HUid]), SIMPLIFY = FALSE);
    ResultMat <- Reduce(f = '+', ResultMatList) + IdMat;
    
    return(ResultMat)
}


## Function that calculate the combination of kernel matrices in the inner layer
CalcInnerKernelComb <- function(lambdajVec, innerKernelList)
{
    n <- dim(innerKernelList[[1]])[1];  # number of individuals;
    J <- length(lambdajVec);   # number of inner layer kernel matrices
    
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
CalcDerivLoss <- function
(
    lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3)
{
    y <- trait

    n <- NROW(y)                        # sample size
    L <- nrow(lambdaliMat)              # number of base kernels
    m <- ncol(lambdaliMat)              # number of hidden units (U)
    J <- NROW(lambdajVec)               # number of inner kernels

    ## Initialization variables
    Amat <- DevAPhi <- matrix(0, nrow = n, ncol = n); # results for matrix A

    ## results for derivative of A and A^2 wrt lambda's
    DevALambdajList <- replicate(J, matrix(.0, n, n), simplify=F)
    DevALambdaliList <- replicate(L * m, matrix(.0, n, n), simplify=F)
    UMat <- matrix(0, nrow = n, ncol = m);
    
    ## Choose the first L and first J kernels in baseKernelList and innerKernelList respectively
    baseKernelList <- baseKernelList[1:L];
    innerKernelList <- list();
    innerKernelName <- innerKernelName[1:J];
    
    ## Calculate the covariance matrices for the m hidden units (DOES NOT DEPEND ON U MATRIX!)
    bkn <- do.call(cbind, lapply(c(baseKernelList, list(diag(n))), as.vector))
    dim(bkn) <- c(n, n, (L + 1))
    bknCmb <- KW(bkn, rbind(exp(lambdaliMat), 1))
    ## hidden U's variance covariance
    UCV <- apply(bknCmb, 3, FastInverseMatrix) 
    dim(UCV) <- dim(bknCmb)
    
    ## Do sampling here to calculate expectations
    dA.inr <- array(0, c(n, n, length(lambdajVec)))
    dA.bas <- array(0, c(n, n, dim(lambdaliMat)))
    dA.phi <- matrix(0, n, n)
    for(s in 1:nSamp)
    {
        ## Sample U from the multivariate normal distribution
        UMat <- apply(UCV, 3L, function(v) mvrnorm(1, rep(0, n), v))
        trPhi <- sum(UKV(UMat, UCV))
        
        ## Obtain the inner layer kernel matrices
        ## for (j in 1:J)
        ##     innerKernelList[[j]] <- findKernel(innerKernelName[j], geno = UMat);
        innerKernelList <- lapply(innerKernelName, findKernel, geno=UMat)
        ikn <- do.call(cbind, lapply(c(innerKernelList, list(diag(n))), as.vector))
        dim(ikn) <- c(n, n, (J + 1))

        ## inversed re-combination of inner kernels serves as predicted VCV of Y
        YCV <- FastInverseMatrix(KW(ikn, rbind(as.matrix(exp(lambdajVec)), 1)))

        ## For matrix A
        Amat <- (1/nSamp) * (Amat + YCV);
        
        ## For derivatives of A with respect to lambda_j
        tmp <- YCV %*% YCV
        for (j in 1:J)
        {
            DevALambdajList[[j]] <- (1/nSamp) * (
                DevALambdajList[[j]] - exp(lambdajVec[j]) * innerKernelList[[j]] %*% tmp)
        }
        for(j in 1:J)
        {
            dA.inr[, , j] <- (1/nSamp) * (dA.inr[, , j] - exp(lambdajVec[j]) * ikn[, , j] %*% tmp)
        }
        
        ## For derivatives of A with respect to lambda_{li}
        for(k in 1:(L*m))
        {
            l <- getIndex(k, L, m)[1];
            i <- getIndex(k, L, m)[2];
            trMat <- UCV[, , i] %*% bkn[, , l] %*% UCV[, , i] %*%
                (bknCmb[, , i] - exp(-phi) * tcrossprod(UMat[,i]));
            trVal <- tr(trMat);
            tmp <- exp(lambdaliMat[l, i]) / 2 * trVal
            DevALambdaliList[[k]] <- (1/nSamp) * (DevALambdaliList[[k]] - tmp * YCV)
        }
        ret <- mapply(function(l, i)
        {
            trMat <- UCV[, , i] %*% bkn[, , l] %*% UCV[, , i] %*%
                (bknCmb[, , i] - exp(-phi) * tcrossprod(UMat[,i]))
            exp(lambdaliMat[l, i]) / 2 * tr(trMat)
        }, rep(1:L, t=m), rep(1:m, e=L))
        dim(ret) <- c(1, 1, L, m)
        dA.bas <- (1/nSamp) * (dA.bas - kronecker(ret, YCV))

        ## For derivative of phi
        DevAPhi <- (1/nSamp) * (DevAPhi - (1/2) * (m * n - exp(-phi) * trPhi) * YCV);
        dA.phi <- (1/nSamp) * (dA.phi - (1/2) * (m * n - exp(-phi) * trPhi) * YCV)
    }

    derivLoss <- sapply(c(DevALambdajList, DevALambdaliList, list(phi=DevAPhi)), function(d)
    {
        t(y) %*% (d %*% Amat + Amat %*% d) %*% y
    })

    dvt <- list()
    dvt$inr <- apply(dA.inr, 3, function(d)
    {
        t(y) %*% (d %*% Amat + Amat %*% d) %*% y
    })
    dvt$bas <- apply(dA.bas, 3:4, function(d)
    {
        t(y) %*% (d %*% Amat + Amat %*% d) %*% y
    })
    dvt$phi <- drop(t(y) %*% (dA.phi %*% Amat + Amat %*% dA.phi) %*% y)

    list(Amat=Amat, dvt=dvt, derivLoss=derivLoss)
}

## Function for gradient descent
## niter is the maximum iterations for gradient descent
## lr is a fixed learning rate used in the option specified and the first iteration
## in the option of B-Bmethod
GradDesc <- function
(
    lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName,
    trait, nSamp = 1e3, niter = 100, lr = 0.001, tol = 1e-5)
{
    derivList <- CalcDerivLoss(
        lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp);

    y <- trait;
    n <- length(trait)                  # sample size
    L <- NROW(lambdaliMat)              # number of base kernels
    m <- NCOL(lambdaliMat)              # number of hidden units (U)
    J <- NROW(lambdajVec)               # number of inner kernels

    ## Initial Gradients
    lossDeriv <- derivList$derivLoss;
    ldv <- list()
    ldv[[2]] <- derivList$dvt
  
    ## Initial parameters in a vector
    paraVec <- c(c(lambdajVec), c(t(lambdaliMat)), phi)
    lpr <- list()
    lpr[[1]] <- list(inr=lambdajVec, bas=lambdaliMat, phi=phi)

    ## Initial Learning Rate
    lr <- getLearningRate(lr=lr, type = "Specified")
    
    ## Matrix used to store the iteration result
    paraMat <- matrix(0, nrow=niter, ncol=length(paraVec))
    paraMat[1,] <- paraVec;
    for(i in 2:niter)
    {
        cat("Iteration", i, "\n");
        if(i == 2)
        {
            lr <- getLearningRate(lr = lr, type = "Specified");
            paraMat[i, ] <- paraMat[i-1, ] - lr * lossDeriv;
            lpr[[i]] <- mapply(function(p, d) p - lr * d, lpr[[1]], ldv[[i]])
            devFxn <- lossDeriv;
        }
        else
        {
            xn <- paraMat[i-1, ]; xn1 <- paraMat[i-2, ];
            devFxn1 <- devFxn;
            
            derivListIter <- CalcDerivLoss(
                lambdajVec = lpr[[i-1]]$inr,
                lambdaliMat = lpr[[i-1]]$bas,
                phi = lpr[[i-1]]$phi,
                baseKernelList = baseKernelList, innerKernelName = innerKernelName,
                trait = trait, nSamp = nSamp);

            lossDerivInter <- derivListIter$derivLoss;
            devFxn <- lossDerivInter;
      
            lr <- getLearningRate(xn, xn1, devFxn, devFxn1, type = "B-BMethod")
            ldv[[i]] <- derivListIter$dvt
            lr1 <- getLearningRate(
                unlist(lpr[[i-1]]), unlist(lpr[[i-2]]),
                unlist(ldv[[i]]), unlist(ldv[[i-1]]), type = "B-BMethod")
            paraMat[i,] <- paraMat[i-1,] - lr * lossDerivInter;

            lpr[[i]] <- mapply(function(p, d) p - lr * d, lpr[[i-1]], ldv[[i]])
        }
        
        cat("LR = ", lr, "\n");
        cat("lmd.j=[", paraMat[i, 1:J], "], lmd.li=[", paraMat[i, (J+1):(J+L*m)],
            "], phi=", paraMat[i, length(paraVec)], "\n");

        if(sum(abs(paraMat[i,] - paraMat[i-1,])) < tol)
        {break;}
    }

    print(all.equal(lpr[[i]]$inr, paraMat[i, 1:J]))
    print(all.equal(lpr[[i]]$bas, matrix(paraMat[i, (J+1):(J+L*m)], nrow = L, byrow = T)))
    print(all.equal(lpr[[i]]$phi, paraMat[i, length(paraVec)]))
    gdl <- list(lambdajVec = lpr[[i]]$inr, lambdaliMat = lpr[[i]]$bas, phi = lpr[[i]]$phi)
    returnList <- c(derivList, gdl)
    return(returnList)
}


# Function for prediction error
CalcPredErr <- function
(
    lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName,
    trait, nSamp = 1e3, niter = 100, lr = 0.001, tol = 1e-5)
{
    y <- trait;
    gradDescList <- GradDesc(
        lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp, niter, lr, tol);
    
    lambdajVecOpt <- gradDescList$lambdajVec;
    lambdaliMatOpt <- gradDescList$lambdaliMat;
    phiOpt <- gradDescList$phi;
    
    optValList <- CalcDerivLoss(
        lambdajVec = lambdajVecOpt, lambdaliMat = lambdaliMatOpt, phi = phiOpt,
        baseKernelList = baseKernelList, innerKernelName = innerKernelName,
        trait = trait, nSamp = nSamp);
    
    ## Calculate prediction error
    optAMat <- optValList$Amat;
    predErr <- t(y) %*% optAMat %*% y;
    
    optimList <- list(
        lambdajVec = lambdajVecOpt, lambdaliMat = lambdaliMatOpt, phi = phiOpt, predErr = predErr);
    returnList <- c(optValList, optimList);
    return(returnList)
}
