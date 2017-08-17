## Function for getting weights of rare variants
getWeights <- function(tag, geno, maf)
{
    if(tag==1)
    {
        wt<-1;
    }
    if(tag == 2)
    {
        wt<-dbeta(maf,1,25)^2 ;
    }
    if(tag == 3)
    {
        wt<-1 / (maf * (1 - maf));
    } 
    if(tag == 0 || tag == 4)
    {
        wt<--log10(maf);
    }
    return(wt);
}

## get best predictor under LMM
## D is the covariance matrix for random effect;
## Z is the design matrix for random effect;
## R is the covariance matrix for error;
getBPLMM <- function(y, D, Z, mu.y = 0, mu.rand = 0, R = 0)
{
    . <- "get best predictor under LMM \
          D is the covariance matrix for random effect \
          Z is the design matrix for random effect \
          R is the covariance matrix for error"
    C <- Z %*% D;
    V <- Z %*% D %*% t(Z) + R;
    
    BP <- mu.rand + C %*% solve(V) %*% (y - mu.y);
    
    return(BP)
}

## Function for getting covariance matrix (Kernels)
## fromBaseKernel is a string vector indicating the kernel to choose from the base kenels.
getKernelList <- function(fromBaseKernel , UserDef = NULL, geno, maf = NULL, weights = 1)
{
    kernelList <- list();
    baseKernelName <- c('CAR', 'identity', 'product');
    
    TrueIndx <- baseKernelName %in% fromBaseKernel;
    
    for(i in 1:sum(TrueIndx))
    {
        if(isTRUE(TrueIndx[i]))
        {
            kernelName <- baseKernelName[match(fromBaseKernel[i], baseKernelName)];
            kernelList[[kernelName]] <- findKernel(kernelName, geno, maf, weights);
        }
    }
    
    ## append user defined kernels into the kernelList
    if(! is.null(UserDef))
    {
        newLength <- length(UserDef);
        for (i in 1:newLength)
        {
            kernelList[[i + length(kernelList)]] <- UserDef[[i]];
        }
    }
    
    return(kernelList)
}


## Main function for simulation
## tag is to indicate the weight used for rare variants;
## func.frq is the frequency of functional variants;
## nS is the number of simulations to be run;
## nN is the total number of individuals;
getPred <- function
(
    tag, func.frq, nS, nSeq, nN, info, ped, # variables to generate SNP data
    traitMu, traitSigma, traitFun, 
    sigmaR, phi, order,                     # variables to generate trait
    fromBaseKernel, UserDef = NULL,         # variables to generate base kernel matrix list
    lambdajVec, lambdaliMat, varPhi,
    innerKernel, 
    nSamp = 1e3, niter = 100, tol = 1e-5)   # variables to assign initial values to the parameters
{
    predErr <- rep(0, nS);
    predErr.BLUPK <- rep(0, nS);
    predErr.BLUPI <- rep(0, nS);
    
    IdMat <- diag(rep(1, nN));
    
    set.seed(525)
    seg.pos <- runif(nS, 0, 1e6-nSeq);
    n.rare <- NULL;
    
    ## pb <- txtProgressBar(min = 0, max = nS, style = 3); cat("\n");
    for(i in 1:nS)
    {
        idx <- (info$pos > seg.pos[i] & info$pos < seg.pos[i] + nSeq); # get a genomic region
        smp <- sample(1:nrow(ped), nN);     # sample nN individuals
        
        geno <- ped[smp, idx]; # obtain their genetic information within the region
        maf <- info$maf[idx];
        
        vrt <- apply(geno, 2, variant);
        geno <- as.matrix(geno[,vrt]);
        maf <- maf[vrt];
        n.rare <- c(n.rare, ncol(geno));
        
        
        idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq);
        geno.true<-as.matrix(geno[,idx.true]);
        maf.true<-maf[idx.true];
        
        weights <- getWeights(tag = tag, geno = geno.true, maf = maf.true);
        
        ## Generate base kernel matrix list
        KList <- getKernelList(
            fromBaseKernel = fromBaseKernel, UserDef = UserDef,
            geno = geno.true, maf = maf.true, weights = weights);
        
        ## Generate trait list
        traitList <- traitSim(
            mu = traitMu, Sigma = traitSigma,
            sigmaR = sigmaR, phi = phi, order = order, fun = traitFun);
        trait <- traitList$trait;
        
        ## Prediction error calculated here
        ctx <- list(
            y=trait,
            nnt=list(
                bas=list(knl=KList, out=3),
                inr=list(knl=innerKernel)))
        err <- CalcPredErr(ctx, nSamp = nSamp, niter = niter, tol = tol)
        
        predVal.BLUPK <- getBPLMM(y = trait, D = sigmaR * traitSigma, Z = IdMat, R = phi * IdMat);
        predVal.BLUPI <- getBPLMM(y = trait, D = sigmaR * IdMat, Z = IdMat, R = phi * IdMat);
        
        predErr[i] <- err
        predErr.BLUPK[i] <- t(trait - predVal.BLUPK) %*% (trait - predVal.BLUPK);
        predErr.BLUPI[i] <- t(trait - predVal.BLUPI) %*% (trait - predVal.BLUPI);
        
        ## setTxtProgressBar(pb, i);
        cat("\n");
    }
    ## close(pb)
    set.seed(NULL);
    
    result <- matrix(
        c(mean(predErr), sd(predErr),
          mean(predErr.BLUPK), sd(predErr.BLUPK),
          mean(predErr.BLUPI), sd(predErr.BLUPI)), nrow = 1);
    colnames(result) <- c(
        "Pred.Err.Mean", "Pred.Err.SD",
        "Pred.Err.Ker.Mean", "Pred.Err.Ker.SD",
        "Pred.Err.Id.Mean", "Pred.Err.Id.SD");
    return(result)
}
