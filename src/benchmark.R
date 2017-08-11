benchmark <- function()
{
    bkl = replicate(100, matrix(rnorm(100*100), 100, 100), simplify=F)
    lmd = rnorm(100)

    library(abind)
    arr = unname(do.call(abind, list(bkl, along=0)))
    ftb <- arr
    dim(ftb) <- c(100, 100*100)

    time1 <- proc.time()
    for(i in 1:2e5)
    {
        ##    user  system elapsed 
        ## 653.136   0.284 654.133 
        abc <- mapply('*', bkl, lmd, SIMPLIFY = FALSE)
        abc <- Reduce('+', abc)
    }
    time2 <- proc.time()
    for(i in 1:2e5)
    {
        ##    user  system  elapsed 
        ## 413.848  94.284  508.119 
        def <- lmd * arr
        def <- colSums(def)
    }
    time3 <- proc.time()
    for(i in 1:2e5)
    {
        ##    user    system  elapsed 
        ## 523.684  1006.736  195.166 
        ghi <- crossprod(lmd, ftb)
        dim(ghi) <- c(100, 100)
    }
    time4 <- proc.time()

    print(all.equal(abc, def))
    print(all.equal(abc, ghi))
    print(time2 - time1)
    print(time3 - time2)
    print(time4 - time3)
}

kW <- function(k, w)
{
    ## an analogy of a left vector times a right matrix, where the N-N-K kernel
    ## array is seen as a row vector of P kernels of size N, which is to be
    ## mapped to another row vector of Q kernels of the same size via Q sets of
    ## weights in the P-Q transformation matrix {w}.
    ## inputs:
    ## k: a 'row' of K input kernels of size N, organized into a N-N-K array;
    ## w: P sets of recombination weights put in a K-P matrix;
    ##
    ## return: P new kernels of size N, arranged in a N-N-P array, where the i
    ## th kernel R[, , i] is the sum of K input kernels weighted by w[, i].
    dk <- dim(k)
    dw <- dim(w)
    stopifnot(tail(dk, 1) == head(dw, 1))

    ## flatten the N-N-K kernal array into N^2-K
    dim(k) <- c(prod(head(dk, -1)), tail(dk, 1))

    ## get P recombinations
    k <- k %*% w                        # N^2-K \times K-P -> N^2-P

    ## restore the N-N kernel dimenstions
    dim(k) <- c(head(dk, -1), tail(dw, 1))

    ## return
    k
}

WKW <- function(k, w)
{
    ## calculate the K quadratic sum w_i' k_i w_i for all i=1...K kernels.
    ##
    ## inputs:
    ## k: K kernels of size N, organized in a N-N-K array;
    ## w: K sets of N-weights put in a N-K matrix;
    ##
    ## return: K summation in a vector, where the i th. entry is the quadratic
    ## sum of the k th. kernel k_i weighted by w_i:
    ## r_i = w_i' k_i w_i
    K = dim(k)[3]
    r = double(K)
    for(i in 1:K) r[i] <- crossprod(w[, i], k[, , i] %*% w[, i])
    r
}

bmk1 <- function(P=4, N=1000)
{
    ## 3-axis array of P kernels (thay are symmetric)
    arr <- rnorm(N^2 * P); dim(arr) <- c(N, N, P)
    knl <- list()
    for(i in 1:P)
    {
        arr[, , i] <- crossprod(arr[, , i])
        knl <- c(knl, list(arr[, , i]))
    }

    ## weight matrix
    lmd <- matrix(rnorm(N * P), N, P)
    print('get ready')

    ## method 1: loops
    t0 <- proc.time()
    rt1 <- sum(sapply(1:P, function(i) crossprod(lmd[, i], arr[, , i] %*% lmd[, i])))
    t1 <- proc.time()
    print(t1 - t0)

    t0 <- proc.time()
    rt3 <- sum(WKW(arr, lmd))
    t1 <- proc.time()
    print(t1 - t0)
    
    ## method 2: lists and traces
    t1 <- proc.time()
    trs <- 0.0
    for(i in 1:P) trs <- trs + sum(diag(arr[, , i] %*% tcrossprod(lmd[, i])));
    rt2 <- trs
    t2 <- proc.time()
    print(t2 - t1)

    print(all.equal(rt1, rt2))
}

bmk3 <- function(P=20, Q=40, N=2000)
{
    ## 3-axis array format
    arr <- rnorm(N^2 * P)
    dim(arr) <- c(N, N, P)

    ## list of matrix format
    knl <- lapply(1:P, function(i) arr[, , i])

    ## weights
    lmd <- matrix(rnorm(P * Q), P, Q)

    ## list routine
    cmb <- list()
    t1 <- proc.time()
    for (i in 1:Q)
    {
        rst <- mapply(`*`, knl, lmd[, i], SIMPLIFY = FALSE);
        cmb[[i]] <- Reduce(`+`, rst)
    }
    t2 <- proc.time()
    arr <- kW(arr, lmd)
    t3 <- proc.time()
    for(i in 1:Q)
    {
        print(all.equal(cmb[[i]], arr[, , i]))
    }
    print(t2-t1)
    print(t3-t2)
}
