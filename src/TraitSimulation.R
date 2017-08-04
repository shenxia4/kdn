# Function Used to Simulate Different Traits Based on Different Setup
# loading necessary libraries
# library(MASS)

variant<-function(x)
{
  length(unique(x)) > 1;
}


traitSim <- function(mu, Sigma, sigmaR = 1, phi = 1, order = 1, fun = c("identity", "polynomial"))
{
  nN <- length(mu);   # number of individuals
  IdMat <- diag(rep(1, nN));
  
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
  
  return(list(trait = y, sigmaR = sigmaR, phi = phi))
  
}