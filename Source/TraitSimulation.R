# Function Used to Simulate Different Traits Based on Different Setup
# loading necessary libraries
# library(MASS)

variant<-function(x)
{
  length(unique(x)) > 1;
}


traitSim <- function(mu, SigmaName, geno, maf, weights, sigmaR = 1, phi = 1, order = 1, fun = c("identity", "polynomial", "sin"))
{
  nN <- length(mu);   # number of individuals
  IdMat <- diag(rep(1, nN));
  
  Sigma <- findKernel(SigmaName, geno = geno, maf = maf, weights = weights)
  
  #set.seed(525)
  ranEff <- mvrnorm(n = 1, mu = mu, Sigma = sigmaR * Sigma);
  
  if(fun == "identity")
  {
    y <- ranEff + sqrt(phi) * rnorm(nN);
  }
  if(fun == "polynomial")
  {
    ranEff <- ranEff ^ order;
    y <- ranEff + sqrt(phi) * rnorm(nN);
  }
  if(fun == "sin")
  {
    ranEff <- sin(ranEff);
    y <- ranEff + sqrt(phi) * rnorm(nN);
  }
  
  return(list(trait = y, sigmaR = sigmaR, phi = phi))
  
}