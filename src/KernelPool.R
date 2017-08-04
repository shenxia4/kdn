## main function
findKernel <- function(kernelName, geno, maf = NULL, weights = 1)
{
  switch(kernelName,
    CAR = {getCAR(geno = geno, maf = maf, weights = weights)},
    identity = {getIdentity(geno = geno);},
    product = {getProduct(geno = geno);},
    stop("Invalied Kernel Name!")
  )
}


## Folllowings are the definition of base kernels
getIBS <- function(geno,weights = 1)
{
  n <- nrow(geno);m <- ncol(geno);
  gtemp1 <- geno;gtemp2 <- geno;
  gtemp1[geno == 2] <- 1;gtemp2[geno == 1] <- 0;gtemp2[geno == 2]<- 1;
  gtemp <- cbind(gtemp1,gtemp2);
  Inner <- gtemp %*% diag(weights,nrow = 2*m,ncol = 2*m) %*% t(gtemp);
  X2 <- matrix(diag(Inner),nrow = n,ncol = n);
  Y2 <- matrix(diag(Inner),nrow = n,ncol = n,byrow = T);
  Bound <- sum(matrix(weights,nrow = 1,ncol = m) * 2);
  Dis <- (X2 + Y2 - 2 * Inner);
  IBS <- Bound - Dis;
  IBS;
}


getCAR <- function(geno, maf, weights)
{
  ## CAR Kernel
  S<-getIBS(geno,weights)/(2*sum(weights));
  if(weights==1){S<-getIBS(geno,weights)/(2*ncol(geno))};
  diag(S)<-0;
  diag<-rowSums(S);
  D<-diag(diag);
  
  gamma <- mean(cor(geno));
  
  Va <- D - gamma * S;
  VaL <- chol(Va);
  K <- chol2inv(VaL);
  
  return(K)
}

getIdentity <- function(geno)
{
  ## identity Kernel
  n <- nrow(geno);
  K <- diag(rep(1, n));
  return(K)
}

getProduct <- function(geno)
{
  K <- geno %*% t(geno);
  return(K)
}


