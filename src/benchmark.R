## bm2hol <- function()
## {
##     bkl = replicate(100, matrix(rnorm(100*100), 100, 100), simplify=F)
##     lmd = rnorm(100)

##     library(abind)
##     arr = unname(do.call(abind, list(bkl, along=0)))
##     ftb <- arr
##     dim(ftb) <- c(100, 100*100)

##     time1 <- proc.time()
##     for(i in 1:2e5)
##     {
##         ##    user  system elapsed 
##         ## 653.136   0.284 654.133 
##         abc <- mapply('*', bkl, lmd, SIMPLIFY = FALSE)
##         abc <- Reduce('+', abc)
##     }
##     time2 <- proc.time()
##     for(i in 1:2e5)
##     {
##         ##    user  system  elapsed 
##         ## 413.848  94.284  508.119 
##         def <- lmd * arr
##         def <- colSums(def)
##     }
##     time3 <- proc.time()
##     for(i in 1:2e5)
##     {
##         ##    user    system  elapsed 
##         ## 523.684  1006.736  195.166 
##         ghi <- crossprod(lmd, ftb)
##         dim(ghi) <- c(100, 100)
##     }
##     time4 <- proc.time()

##     print(all.equal(abc, def))
##     print(all.equal(abc, ghi))
##     print(time2 - time1)
##     print(time3 - time2)
##     print(time4 - time3)
## }

KW <- function(k, w, drop=TRUE)
{
    ## an analogy of transforming a K-vector {k} by a K-P weight matrix {w}, such
    ## that r_j = k'w[,j], for j=1,...,P, where a list of K kernels are seen as K
    ## entries of the K-vector {k}.
    ## 
    ## inputs:
    ## {k}: K input kernels of size N, organized into a N-N-K array;
    ## {w}: P sets of recombination weights put in a K-P matrix;
    ##
    ## return: P new kernels of size N, arranged in a N-N-P array, where the j th
    ## r[, , j] is the sum of K input kernels weighted by w[, j].
    dk <- dim(k)
    dw <- dim(w)
    stopifnot(tail(dk, 1) == head(dw, 1))

    ## flatten the N-N-K kernal array into N^2-K
    dim(k) <- c(prod(head(dk, -1)), tail(dk, 1))

    ## get P recombinations
    k <- k %*% w                        # N^2-K \times K-P -> N^2-P

    ## restore the N-N kernel dimenstions
    dim(k) <- c(head(dk, -1), tail(dw, 1))

    ## drop/squeeze axis of 1 dimension?
    if(drop)
        k <- drop(k)
    k
}

bm1 <- function(n=50, s=5000)
{
    n^2 %>% rnorm %>% matrix(n) %>% crossprod -> x
    u <- chol(x)
    v <- solve(x, diag(n))
    l <- t(u)
    
    ## facts:
    ## l == t(u) and t(l) == u
    ## l %*% u == l %*% t(l) == t(u) %*% u == x
    ## benchmark(
        ## iu <- backsolve(u, I, n, T),    # ***
        ## il <- backsolve(l, I, n, F),    # **
        ## ic <- chol2inv(u),              # *
        ## replications=5)
    ## facts: iu %*% t(iu) == t(il) %*% il) == inv(x)

    (s*n) %>% rnorm %>% matrix(s, n) -> s
    ## benchmark(
    ##     sc <- a %*% chol(ic),
    ##     su <- tcrossprod(a, iu),
    ##     replications=5)

    ## benchmark(
    ##     sc <- s %*% chol(chol2inv(u)),
    ##     su <- tcrossprod(s, backsolve(u, diag(n))),
    ##     replications=5) %>% print

    benchmark(
        sc <- s %*% chol(chol2inv(u)),
        su <- tcrossprod(s, backsolve(u, diag(n))),
        replications=2) %>% print
    list(v=v, sc=sc, su=su)
}
