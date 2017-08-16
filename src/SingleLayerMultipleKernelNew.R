## source("HelperFunctions.R")
## source("KernelPool.R")

## Function used to calculate necessary derivatives for gradient descent algorithm
## lambdajVec is the initial vector for the variance components for final prediction;
## lambdaliMat is an L by m matrix containing the initial values for the varaince component for the
## single layer;
## L is the number of base kernel matrices;
## J is the number of inner layer kernel matrices;
## m is the number of hidden units;
## trait is the phenotype vector;
## baseKernelList is a matrix list containing the L base kernel matrices;
## innerKernelName list the name of all the kernels to be used in the inner layer;
## baseKernelList can have more than L elements and innerKernelName can have more than J elements.
## If so, the first L and J elements will be used;
CalcDerivLoss <- function
(
    lambdajVec, lambdaliMat, phi, baseKernelList, innerKernelName, trait, nSamp = 1e3)
{
    y <- trait

    N <- NROW(y)                        # sample size
    L <- nrow(lambdaliMat)              # number of base kernels
    m <- ncol(lambdaliMat)              # number of hidden units (U)
    J <- NROW(lambdajVec)               # number of inner kernels

    ## Initialization variables
    Amat <- matrix(0, N, N);            # A matrix
    UMat <- matrix(0, N, m);            # U matrix (hidden units)

    ## Choose the first L and first J kernels in baseKernelList and innerKernelList respectively
    baseKernelList <- baseKernelList[1:L];
    innerKernelList <- list();
    innerKernelName <- innerKernelName[1:J];

    ## Calculate the covariance matrices for the m hidden units (DOES NOT DEPEND ON U MATRIX!)
    bkn <- do.call(cbind, lapply(c(baseKernelList, list(diag(N))), as.vector))
    dim(bkn) <- c(N, N, (L + 1))
    bknCmb <- KW(bkn, rbind(exp(lambdaliMat), 1))
    ## hidden U's variance covariance
    UCV <- apply(bknCmb, 3, FastInverseMatrix)
    dim(UCV) <- dim(bknCmb)

    ## gradient of A wrt weights
    dA.inr <- array(0, c(N, N, length(lambdajVec)))
    dA.bas <- array(0, c(N, N, dim(lambdaliMat)))
    dA.phi <- matrix(0, N, N)

    ## y' \PDV{A}{\lambda_j }, for all j,   -- gradient of A wrt. inner weights times y
    dA.inr.y <- array(0, c(N, J))

    ## \PDV{A}{\lambda_li}y, for all (l, i) -- gradient of A wrt. basic weights times y
    dA.bas.y <- array(0, c(N, L, m))

    ## y' \PDV{A}{\phi}, y' time the derivative of A wrt. phi
    dA.phi.y <- 0
    
    for(s in 1:nSamp)
    {
        ## Sample U from the multivariate normal distribution
        UMat <- apply(UCV, 3L, function(v) mvrnorm(1, rep(0, N), v))
        
        ## Obtain the inner layer kernel matrices
        innerKernelList <- lapply(innerKernelName, findKernel, geno=UMat)
        ikn <- do.call(cbind, lapply(c(innerKernelList, list(diag(N))), as.vector))
        dim(ikn) <- c(N, N, (J + 1))
        
        ## inversed inner kernels combination serves as predicted VCV of Y,
        ## which is also a sample of matrix A: YVC = A[s]
        YCV <- FastInverseMatrix(KW(ikn, rbind(as.matrix(exp(lambdajVec)), 1)))

        ## accmulate matrix A
        Amat <- Amat + YCV
        
        ## For derivatives of A with respect to inner weights lambda_j
        .yy <- YCV %*% YCV
        dA.inr <- dA.inr + vapply(1:J, function(j)
        {
            -exp(lambdajVec[j]) * ikn[, , j] %*% .yy
        }, YCV)
        dA.inr.y <- dA.inr.y + vapply(1:J, function(j)
        {
            ## \dev(A, lmd_j) y = A K_j A y = y' A K_j A for all j.
            -exp(lambdajVec[j]) * crossprod(y, YCV) %*% ikn[, , j] %*% YCV
        }, y)

        ## For derivatives of A with respect to basic weights lambda_{li}
        ret <- mapply(function(l, i)
        {
            . <- UCV[, , i] %*% UMat[, i]
            . <- - exp(-phi) * crossprod(., bkn[, , l] %*% .)
            . <- sum(UCV[, , i] * bkn[, , l]) + .
            -.5 * exp(lambdaliMat[l, i]) * .
        }, rep(1:L, t=m), rep(1:m, e=L))
        dim(ret) <- c(1, 1, L, m)
        dA.bas <- dA.bas + kronecker(ret, YCV)

        ## \dev(A, lmd_li)y = A [..] y = y' A [..] for all (l, i)
        dim(ret) <- c(1, L, m)
        dA.bas.y <- dA.bas.y + kronecker(ret, drop(crossprod(y, YCV)))
        
        ## For derivative of phi
        . <- sum(sapply(1:m, function(i) crossprod(UMat[, i], UCV[, , i] %*% UMat[, i])))
        dA.phi <- dA.phi - .5 * YCV * (m * N - exp(-phi) * .)
        ## \dev(A, phi) y = A . y = (y' A)' . = A y .
        dA.phi.y <- dA.phi.y + YCV %*% y * -.5 * (m * N - exp(-phi) * .)
    }
    Amat <- Amat / nSamp
    dA.inr <- dA.inr / nSamp
    dA.bas <- dA.bas / nSamp
    dA.phi <- dA.phi / nSamp

    t1 <- proc.time()
    dA.inr.y <- dA.inr.y / nSamp
    dA.inr.2 <- drop(crossprod(y, Amat %*% dA.inr.y)) * 2

    dA.bas.y <- dA.bas.y / nSamp
    dim(dA.bas.y) <- c(N, L * m)       # flatten
    dA.bas.2 <- crossprod(y, Amat %*% dA.bas.y) * 2
    dim(dA.bas.2) <- c(L, m)           # reshape

    dA.phi.y <- dA.phi.y / nSamp
    dA.phi.y <- drop(crossprod(y, Amat) %*% dA.phi.y) * 2
    t2 <- proc.time()
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
    t3 <- proc.time()
    print(t2-t1)
    print(t3-t2)
    
    list(Amat=Amat, dvt=dvt)
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
    N <- length(trait)                  # sample size
    L <- NROW(lambdaliMat)              # number of base kernels
    m <- NCOL(lambdaliMat)              # number of hidden units (U)
    J <- NROW(lambdajVec)               # number of inner kernels

    ## Initial Gradients
    ldv <- list()
    ldv[[2]] <- derivList$dvt
  
    ## Initial parameters in a vector
    lpr <- list()
    lpr[[1]] <- list(inr=lambdajVec, bas=lambdaliMat, phi=phi)

    ## Initial Learning Rate
    lr <- getLearningRate(lr=lr, type = "Specified")
    
    ## Matrix used to store the iteration result
    for(i in 2:niter)
    {
        cat("Iteration", i, "\n");
        if(i == 2)
        {
            lr <- getLearningRate(lr = lr, type = "Specified");
            lpr[[i]] <- mapply(function(p, d) p - lr * d, lpr[[1]], ldv[[i]])
        }
        else
        {
            derivListIter <- CalcDerivLoss(
                lambdajVec = lpr[[i-1]]$inr,
                lambdaliMat = lpr[[i-1]]$bas,
                phi = lpr[[i-1]]$phi,
                baseKernelList = baseKernelList, innerKernelName = innerKernelName,
                trait = trait, nSamp = nSamp);
            ldv[[i]] <- derivListIter$dvt
            lr1 <- getLearningRate(
                unlist(lpr[[i-1]]), unlist(lpr[[i-2]]),
                unlist(ldv[[i]]), unlist(ldv[[i-1]]), type = "B-BMethod")
            lpr[[i]] <- mapply(function(p, d) p - lr1 * d, lpr[[i-1]], ldv[[i]])
        }
        
        cat("LR = ", lr, "\n")
        cat("lmd.j=[", lpr[[i]]$inr, "], lmd.li=[", lpr[[i]]$bas, "], phi=", lpr[[i]]$phi, "\n")
        if(sum(abs(unlist(lpr[[i-1]]) - unlist(lpr[[i]]))) < tol)
            break
    }

    list(lambdajVec = lpr[[i]]$inr, lambdaliMat = lpr[[i]]$bas, phi = lpr[[i]]$phi)
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
