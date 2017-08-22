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
## FastInverseMatrix <- function(X)
## {
##     XChol <- chol(X)
##     XInv <- chol2inv(XChol)
##     return(XInv)
## }

## Trace Function
tr <- function(X)
{
    return(sum(diag(X)))
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
