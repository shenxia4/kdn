## source("HelperFunctions.R")
## source("KernelPool.R")

## Function used to calculate necessary derivatives for gradient descent algorithm
## par: parameters
##    inr: inner weights
##    bas: basic weights
##    phi: phi
## y is the phenotype vector
## knl: kernels
##    bas: basic kernels (matrices)
##    inr: inner kernels (tokens)
CalcDerivLoss <- function(par, knl, y, nSamp = 1e3)
{
    ## dimensions
    N <- NROW(y)                        # sample size
    L <- nrow(par$bas)                  # number of base kernels
    M <- ncol(par$bas)                  # number of hidden units (U)
    J <- NROW(par$inr)                  # number of inner kernels
    K <- NCOL(par$inr)                  # number of output units (Y)

    ## Initialization variables
    Amat <- matrix(0, N, N);            # A matrix
    UMat <- matrix(0, N, M);            # U matrix (hidden units)
    
    ## organize the basic kernels into a 3D array where the last axis indices kernels,
    ## also append an identical kernel
    bkn <- vapply(c(knl$bas[1:L], list(diag(N))), unname, Amat)

    ## get the covariance structure of M hidden units (NOT DEPEND ON U MATRIX!)
    . <- bkn
    dim(.) <- c(N^2, L+1)
    . <- . %*% rbind(exp(par$bas), 1)
    dim(.) <- c(N, N, M)
    UCV <- vapply(1:M, function(i) FastInverseMatrix(.[, , i]), Amat)
    
    ## gradient accumulents
    ## \PDV{A}{\lambda_j }y, for all j,   -- gradient of A wrt. inner weights times y
    dA.inr.y <- array(0, c(N, J))

    ## \PDV{A}{\lambda_li}y, for all (l, i) -- gradient of A wrt. basic weights times y
    dA.bas.y <- array(0, c(N, L, M))

    ## \PDV{A}{\phi}y, -- gradient of A wrt. phi time y
    dA.phi.y <- 0
    
    for(s in 1:nSamp)
    {
        ## Sample U from the multivariate normal distribution
        UMat <- apply(UCV, 3L, function(v) mvrnorm(1, rep(0, N), v))
        
        ## Obtain the inner layer kernel matrices
        ikn <- c(lapply(knl$inr[1:J], findKernel, geno=UMat), list(diag(N)))
        ikn <- vapply(ikn, I, Amat)
        
        ## inversed inner kernels combination serves as predicted VCV of Y,
        ## which is also a sample of matrix A: YVC = A[s]
        . <- ikn
        dim(.) <- c(N^2, J+1)
        . <- . %*% c(exp(par$inr), 1)
        dim(.) <- c(N, N)
        YCV <- FastInverseMatrix(.)
        
        ## accmulate matrix A
        Amat <- Amat + YCV

        ## y'A
        yT.A <- crossprod(y, YCV)
        
        ## For derivatives of A with respect to inner weights lambda_j
        dA.inr.y <- dA.inr.y + vapply(1:J, function(j)
        {
            ## \dev(A, lmd_j) y = A K_j A y = y' A K_j A for all j.
            -exp(par$inr[j]) * yT.A %*% ikn[, , j] %*% YCV
        }, y)

        ## For derivatives of A with respect to basic weights lambda_{li}
        ret <- mapply(function(l, i)
        {
            . <- UCV[, , i] %*% UMat[, i]
            . <- - exp(-par$phi) * crossprod(., bkn[, , l] %*% .)
            . <- sum(UCV[, , i] * bkn[, , l]) + .
            -.5 * exp(par$bas[l, i]) * .
        }, rep(1:L, t=M), rep(1:M, e=L))
        ## \dev(A, lmd_li)y = A [..] y = y' A [..] for all (l, i)
        dim(ret) <- c(1, L, M)
        dA.bas.y <- dA.bas.y + kronecker(ret, drop(yT.A))
        
        ## For derivative of phi
        . <- sum(sapply(1:M, function(i) crossprod(UMat[, i], UCV[, , i] %*% UMat[, i])))
        ## \dev(A, phi) y = A [.] y = (y' A)' . = A y [.]
        dA.phi.y <- dA.phi.y + drop(yT.A) * -.5 * (M * N - exp(-par$phi) * .)
    }

    Amat <- Amat / nSamp
    dA.inr.y <- dA.inr.y / nSamp
    dA.bas.y <- dA.bas.y / nSamp
    dA.phi.y <- dA.phi.y / nSamp

    dim(dA.bas.y) <- c(N, L * M)        # flatten
    dA.bas.2 <- crossprod(y, Amat %*% dA.bas.y) * 2
    dim(dA.bas.2) <- c(L, M)            # reshape

    dvt <- within(list(),
    {
        phi <- drop(crossprod(y, Amat) %*% dA.phi.y) * 2
        bas <- dA.bas.2
        inr <- drop(crossprod(y, Amat %*% dA.inr.y)) * 2
    })
    
    list(Amat=Amat, dvt=dvt)
}

## Function for gradient descent
## ctx     : training contex
##    par  : list of parameters -- inner weights, basic weights, and phi
##    knl  : list of kernels    -- inner kernels, basic kernels
##      y  : response (ie. phenotype)
##
## ---- configurations ----
## max.itr : maximum number of iterations;
## lr      : is the initial learning rate;
## min.err : minimum training error to continue;
## nSamp   : number of sampling to do when evaluating gradient
GradDesc <- function(ctx, nSamp=1e3, max.itr=100, lr=1e-3, tol=1e-5, min.err=1e-3)
{
    for(. in names(ctx)) assign(., ctx[[.]])
    ## Matrix used to store the iteration result
    for(i in 1:max.itr)
    {
        if(i == 1)                      # Initial: 1st iteration
        {
            lpr <- list()               # parameters
            ldv <- list()               # gradients
            lpr[[i]] <- par
        }

        ## calculate gradient and other context variables
        par <- lpr[[i]]
        ret <- with(par, CalcDerivLoss(par, knl, y, nSamp))
        dvt <- ldv[[i]] <- ret$dvt

        ## update parameters
        if(i < 3)                       # Initial learning rate
            lr <- getLearningRate(lr=lr, type = "Specified")
        else                            # dynamic learning rate
        {
            p.i <- unlist(par)
            p.1 <- unlist(lpr[[i-1]])
            d.i <- unlist(dvt)
            d.1 <- unlist(ldv[[i-1]])
            ## lr <- getLearningRate(p.i, p.1, d.i, d.1, type="B-BMethod")
        }
        ## update now
        new <- lpr[[i+1]] <- within(par,
        {
            phi <- phi - lr * ret$dvt$phi
            bas <- bas - lr * ret$dvt$bas
            inr <- inr - lr * ret$dvt$inr
        })

        ## statistic tracks
        ## prediction error by previous parameters (i th.)
        err <- t(y) %*% ret$Amat %*% y

        ## difference
        dff <- sum(abs(unlist(new) - unlist(par)))

        ## maximum gradient
        mdv <- max(abs(unlist(dvt)))

        rpt <- sprintf("%04d: %.1e %.1e %.1e %.1e", i, err, lr, mdv, dff)
        cat(rpt, '\n', sep='')

        if(dff < tol)
            break
        if(err < min.err)
            break
    }
    list(par=new)
}


# Function for prediction error
CalcPredErr <- function
(
    lambdajVec, lambdaliMat, phi, baseKernelList, innerKernel,
    trait, nSamp = 1e3, niter = 100, lr = 0.001, tol = 1e-5, ctx=NULL)
{
    ## initialize gradient decent context
    nnt <- ctx$nnt
    knl <- lapply(nnt, `[[`, 'knl')
    lmd <- lapply(nnt, function(.)
    {
        d.i <- length(.$knl)
        d.o <- .$out
        if(is.null(d.o))
            rep(1, d.i)
        else
            matrix(1:d.o - 1, d.i, d.o, T)
    })
    
    ctx <- within(list(),
    {
        y <- ctx$y
        par <- c(lmd, list(phi=1))
        knl <- knl
    })
    ret <- GradDesc(ctx, nSamp=nSamp, max.itr=niter, lr=lr, min.err=tol);

    ## result
    par <- ret$par
    knl <- list(bas=baseKernelList, inr=innerKernel)
    optValList <- CalcDerivLoss(par=par, knl, y=trait, nSamp=nSamp)
    
    ## Calculate prediction error
    optAMat <- optValList$Amat;
    predErr <- t(trait) %*% optAMat %*% trait;
    
    optimList <- list(lambdajVec=par$inr, lambdaliMat=par$bas, phi=par$phi, predErr=predErr)
    returnList <- c(optValList, optimList);
    return(returnList)
}
