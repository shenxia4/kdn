## Function for Calculating Square Root of a Symmetric Matrix
expMat <- function(X, order = 1/2)
{
    if(!isSymmetric(X))
    {
        print("The input matrix is not symmetric!");  
    }else
    {
        U <- eigen(X)$vectors;
        Lambda <- diag((eigen(X)$values) ^ order);
        R <- U %*% Lambda %*% t(U);
    }
    
    return(R)
}

## Function Used to Fastly Inverse a Large Matrix
FastInverseMatrix <- function(X)
{
    r = try(chol(X))
    if(inherits(r, 'try-error'))
    {
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@6@"]]));##:ess-bp-end:##
        print(r)
        return(r)
    }
    XChol <- r
    XInv <- chol2inv(XChol);
    return(XInv)
}

## Trace Function
tr <- function(X)
{
    return(sum(diag(X)))
}

## Initialize a List of matrices
InitializeList <- function(m = 1, n = 100)
{
    ls <- list();
    for(i in 1:m)
    {
        ls[[i]] <- matrix(0, nrow = n, ncol = n);
    }
    
    return(ls)
}

## Function used to obtain the matrix index from a row-based vectorized sequence
## m, n are the row number and the column number of the matrix respectively
## id is the index from the sequenced vector
getIndex <- function(id, m, n)
{
    if(id %% n == 0){mindx <- id / n;}
    else{mindx <- id %/% n + 1; }
    
    nindx <- (id - (mindx - 1) * n) %% n;
    if(nindx == 0){nindx <- n;}
    
    return(c(mindx, nindx))
}


## Function used to obtain the learning rate
getLearningRate <- function(xn, xn1, devFxn, devFxn1, lr = 0.001, type = c("Specified", "B-BMethod"))
{
    if(type == "Specified")
    {
        lr <- as.numeric(lr);
    }
    
    if(type == "B-BMethod")
    {
        lr <- crossprod(xn - xn1, devFxn - devFxn1) / crossprod(devFxn - devFxn1, devFxn - devFxn1);
        lr <- as.numeric(lr);
    }
    return(lr)
}
