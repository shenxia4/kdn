source("src/HelperFunctions.R")
## source("KernelPool.R")
library('magrittr')

#' Sum of vector elements.
#'
#' \code{.chl.inv.mvn} draws multi-variant normal (mvn) sample whose covariance is
#' the inverse of a symmetric matrix \code{x}.
#'
#' By sharing the upper-triangle, the function avoids double Cholesky decomposition
#' typically needed to speed up the inversion of \code{x}, and the mvn sampling from
#' \code{x}^{-1}
#' 
#' It should roughly halv the total execution time required by the two operations if
#' they were run seperately.
#' 
#' @param x numeric, positive definite matrix. no sanity checking is imposed.
#' @param n numeric, positve scalar, the number of sample mvn samples to be drawn.
#' @param ret.inv, logical, TRUE return inverse of \code{x} as well, otherwise, only
#' perform the mvn sampling
#'
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' n <- 100
#' p <- 500
#' x <- crossprod(matrix(rnorm(p^2), p, p))
#' 
#' y <- .chl.inv.mvn(x, n)
#' 
#' \dontrun{
#' #' sum("a")
#' }
.chl.inv.mvn <- function(x, n=1, ret.inv=FALSE)
{
    ## upper triangle of Cholesky decomposition of x
    u <- chol(x)
    p <- nrow(x)

    ## upper-tri of the the inverse x^{-1}
    v <- backsolve(u, diag(p))
    
    ## mvn sampling from N(I_p, x^{-1})
    mvn <- (n * p) %>% rnorm %>% matrix(n, p) %>% tcrossprod(v)

    ## return inverse of x as well?
    ret <- if(ret.inv) list(mvn=mvn, inv=tcrossprod(v)) else mvn
    ret
}

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
    }
    else
    {
        X <- matrix(rnorm(p * n), n, p)
        ## X <- X %*% (diag(A) %*% t(eS$vectors))
        X <- X %*% (A * t(eS$vectors))
    }
    X
}


#' calculate necessary derivatives for gradient descent algorithm
#' @param par list, parameters to calculate derivatives
#'    inr: numeric inner kernel weights, a vector for now
#'    bas: numeric basic kernel weights, a L by M matrix
#'    phi: numeric the phi scalar
#' @param y numeric, the phenotype vector
#' @param knl list, of numeric matrix, string names, or bi-operand
#' function.
#'    bas: basic kernels (usually matrices)
#'    inr: inner kernels (usually names)
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

    ## parameters
    tao.bas <- exp(par$bas) # weights connecting basic kernels and UCVs
    tao.inr <- exp(par$inr) # weights connecting inner kernels and YCV
    phi <- par$phi          # the 'round' phi
    PHI <- exp(phi)         # the 'stick through circle' phi

    ## organize the basic kernels into a 3D array, with last axis indexing the
    ## kernels, also append an identical kernel at the end
    bkn <- vapply(c(knl$bas[1:L], list(diag(N))), unname, Amat)

    ## help function to alter dimensions, and multiply by matrix
    .dim <- function(x, ...) {dim(x) <- c(...); x}
    .mbm <- .Primitive('%*%')
    ## the covariance of M hidden units U_{1...M} (UCV, NOT DEPEND ON U MATRIX!), is
    ## based on the sum of base kernels {bkn} weighted by {tao.bas}.
    ## 1) the M weighted sum and their Cholesky decomposition
    BMX <- bkn %>% .dim(N^2, L+1) %>% .mbm(rbind(tao.bas, 1)) %>% .dim(N, N, M)
    ## 2) the inversed combination
    IUC <- vapply(1:M, function(i) chol2inv(chol(BMX[, , i])), Amat)
    ## 3) the VCV of M hidden units U_{1...M}
    UCV <- vapply(1:M, function(i) BMX[, , i] * PHI, Amat)
    
    ## gradient accumulents
    ## \PDV{A}{\lambda_j }y, for all j,   -- gradient of A wrt. inner weights times y
    dA.inr.y <- array(0, c(N, J))

    ## \PDV{A}{\lambda_li}y, for all (l, i) -- gradient of A wrt. basic weights times y
    dA.bas.y <- array(0, c(N, L, M))

    ## \PDV{A}{\phi}y, -- gradient of A wrt. phi time y
    dA.phi.y <- 0

    ## sampling of hidden units
    UArr <- vapply(1:M, function(i)
        mvrnorm(nSamp, rep(0, N), UCV[, , i]), matrix(.0, nSamp, N))
    for(s in 1:nSamp)
    {
        ## N by M, the {s} th. sample from mvn(0, UCV[, , i]) for i = 1...M
        UMat <- UArr[s, ,]
        
        ## Obtain the inner layer kernel matrices
        ikn <- c(lapply(knl$inr[1:J], findKernel, geno=UMat), list(diag(N)))
        ikn <- vapply(ikn, I, Amat)
        
        ## inner kernel combination serves as predicted VCV of Y, whose inverse
        ## is also the s th. sample of A matrix: A = IYC = YCV^{-1}
        YCV <- ikn %>% .dim(N^2, J+1) %>% .mbm(c(tao.inr, 1)) %>% .dim(N, N)
        IYC <- chol2inv(chol(YCV))      # inverse of Y covariance
        
        ## accmulate matrix A
        Amat <- Amat + IYC

        ## y'A
        yT.A <- crossprod(y, IYC)
        
        ## For derivatives of A with respect to inner weights lambda_j
        dA.inr.y <- dA.inr.y + vapply(1:J, function(j)
        {
            ## \dev(A, lmd_j) y = A K_j A y = y' A K_j A for all j.
            -exp(par$inr[j]) * yT.A %*% ikn[, , j] %*% IYC
        }, y)

        ## For derivatives of A with respect to basic weights lambda_{li}
        ret <- mapply(function(l, i)
        {
            . <- IUC[, , i] %*% UMat[, i]
            . <- - exp(-phi) * crossprod(., bkn[, , l] %*% .)
            . <- sum(IUC[, , i] * bkn[, , l]) + .
            -.5 * tao.bas[l, i] * .
        }, rep(1:L, t=M), rep(1:M, e=L))
        ## \dev(A, lmd_li)y = A [..] y = y' A [..] for all (l, i)
        dim(ret) <- c(1, L, M)
        dA.bas.y <- dA.bas.y + kronecker(ret, drop(yT.A))
        
        ## For derivative of phi
        . <- sum(sapply(1:M, function(i) crossprod(UMat[, i], IUC[, , i] %*% UMat[, i])))
        ## \dev(A, phi) y = A [.] y = (y' A)' . = A y [.]
        dA.phi.y <- dA.phi.y + drop(yT.A) * -.5 * (M * N - exp(-phi) * .)
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
GradDesc <- function(ctx, max.itr=100, lr=1e-3, min.lr=1e-9, max.lr=1e1, tol=1e-5, min.err=1e-3, ...)
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
        mpv <- max(abs(unlist(par)))          # parameter maximum
        err <- with(ret, y %*% Amat %*% y)    # pred-error

        ## learning rate
        if(i < 1000)                    # Initial learning rate
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
        rpt <- list(ep=i, err=err, lr=lr, mdv=mdv, mpv=mpv, sdf=sdf)
        cat(do.call(sprintf, c("%04d: %.1e %.1e %.1e %.1e %.1e\n", rpt)))

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
CalcPredErr <- function(ctx, niter=100, lr=1e-3, tol=1e-6, ...)
{
    ## initialize gradient decent context
    nnt <- ctx$nnt
    knl <- lapply(nnt, `[[`, 'knl')
    lmd <- lapply(nnt, function(.)
    {
        d.i <- length(.$knl)
        d.o <- .$out
        if(is.null(d.o))
            rep(0, d.i)
        else
            matrix(rnorm(d.i * d.o, 0, .1), d.i, d.o)
    })
    ctx <- within(list(),
    {
        y <- ctx$y
        par <- c(lmd, list(phi=0))
        knl <- knl
    })

    ## perform Gradient Descent
    ret <- GradDesc(ctx, max.itr=niter, lr=lr, min.err=tol, tol=tol, ...)
    
    ## result
    err <- ret$hst[[length(ret$hst)]]$err
    err
}
