## Function Used to Simulate Different Traits Based on Different Setup
## loading necessary libraries
## library(MASS)

variant<-function(x)
{
    length(unique(x)) > 1;
}


traitSim <- function(mu, Sigma, sigmaR = 1, phi = 1, fun=I)
{
    nN <- length(mu);                   # number of individuals
    
    ranEff <- mvrnorm(n = 1, mu = mu, Sigma = sigmaR * Sigma);

    y <- fun(ranEFF) + sqrt(phi) * rnorm(nN);
    y
}
