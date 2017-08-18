## source("HelperFunctions.R")
## source("KernelPool.R")

.mvnm <- function (n=1, mu=0, sigma, tol=1e-06)
{
    p <- nrow(sigma)

    eS <- eigen(sigma, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L]))) 
        stop("'Sigma' is not positive definite")

    A <- sqrt(pmax(ev, 0))              # p-vector
    if (p > n)
    {
        X <- matrix(rnorm(p * n), p, n)
        ## t1 <- system.time(Y <- t(eS$vectors %*% (diag(A) %*% X)))
        ## t2 <- system.time(Y <- t(eS$vectors %*% (A * X)))
        ## t3 <- system.time(eS$vectors %*% diag(A) %*% X)
        ## print(t1)
        ## print(t2)
        ## print(t3)
    }
    else
    {
        X <- matrix(rnorm(p * n), n, p)
        ## X <- X %*% (diag(A) %*% t(eS$vectors))
        X <- X %*% (A * t(eS$vectors))
    }
    X
}

## Function used to calculate necessary derivatives for gradient descent algorithm
## par: parameters
##    inr: inner weights
##    bas: basic weights
##    phi: phi
## y is the phenotype vector
## knl: kernels
##    bas: basic kernels (matrices)
##    inr: inner kernels (tokens)
CalcDerivLoss <- function(par, knl, y, nSamp=1e3, ...)
{
    ## dimensions
    N <- NROW(y)                        # sample size
    L <- nrow(par$bas)                  # number of base kernels
    M <- ncol(par$bas)                  # number of hidden units (U)
    J <- NROW(par$inr)                  # number of inner kernels
    K <- NCOL(par$inr)                  # number of output units (Y)

    ## Initialization variables
    Amat <- matrix(0, N, N)             # A matrix
    UMat <- matrix(0, N, M)             # U matrix (hidden units)
    
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

    print(system.time(UArr <- apply(UCV, 3L, function(v) mvrnorm(nSamp, rep(0, N), v))))
    print(system.time(UArr <- apply(UCV, 3L, function(v) .mvnm(nSamp, rep(0, N), v))))
    dim(UArr) <- c(nSamp, N, M)
    for(s in 1:nSamp)
    {
        ## Sample U from the multivariate normal distribution
        ## UMat <- apply(UCV, 3L, function(v) mvrnorm(1, rep(0, N), v))
        UMat <- UArr[s, ,]
        
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
    
    ## make sure the gradient is in the same order as the parameters 
    dvt <- dvt[names(par)]

    ## return
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
GradDesc <- function(ctx, max.itr=100, lr=1e-3, min.lr=1e-6, max.lr=1e-1, tol=1e-5, min.err=1e-3, ...)
{
    ## contex: hst*, par*, y, knl
    for(. in names(ctx)) assign(., ctx[[.]])

    ## initial epoch, with only the initial parameters
    if(!exists('hst'))
        hst <- list()

    par <- ctx$par
    ## training
    for(i in (1 + length(hst)) : (length(hst) + max.itr))
    {
        ## re-calculate gradient and other contextual variables
        ret <- CalcDerivLoss(par, knl, y, ...)
        dvt <- ret$dvt                        # gradient
        mdv <- max(abs(unlist(dvt)))          # gradient maximum
        err <- with(ret, y %*% Amat %*% y)    # pred-error

        ## learning rate
        if(i < 3)                       # Initial learning rate
            lr <- getLearningRate(lr=lr, type = "Specified")
        else                            # dynamic learning rate
        {
            p.i <- unlist(par)
            p.1 <- unlist(hst[[i-1]]$par)
            d.i <- unlist(dvt)
            d.1 <- unlist(hst[[i-1]]$dvt)
            lr <- getLearningRate(p.i, p.1, d.i, d.1, type="B-BMethod")
            lr <- min(lr, max.lr)
            lr <- max(lr, min.lr)
        }
        ## parameter differences, and its summation
        dff <- lapply(dvt, `*`, lr)
        sdf <- sum(abs(unlist(dff)))

        ## gather report
        rpt <- list(ep=i, err=err, lr=lr, mdv=mdv, sdf=sdf)
        cat(do.call(sprintf, c("%04d: %.1e %.1e %.1e %.1e\n", rpt)))

        ## record history
        hst[[i]] <- c(rpt, list(par=par, dvt=dvt))
        names(hst)[i] <- sprintf("E%04d", i)

        if(sdf < tol || err < min.err)
            break

        ## get new parameters
        par <- mapply(`-`, par, dff)

        ## forget early history of parameters and gradients
        if(i > 2)
        {
            hst[[i-2]]$par <- NULL
            hst[[i-2]]$dvt <- NULL
        }
    }
    ctx$par <- par
    ctx$hst <- hst
    ctx
}


## Function for prediction error
CalcPredErr <- function(ctx, niter=100, lr=0.001, tol=1e-5, ...)
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

    ## perform Gradient Descent
    ret <- GradDesc(ctx, max.itr=niter, lr=lr, min.err=tol, tol=tol, ...)
    
    ## result
    err <- ret$hst[[length(ret$hst)]]$err
    err
}
