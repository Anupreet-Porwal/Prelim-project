library(statmod)
neg.normal.means <- function(y, c,d,burn=1000,nmc=5000){
  
  n <- length(y)
  BetaSave = matrix(0, nmc, n)
  TauSave = rep(0, nmc)
  Sigma2Save = rep(0, nmc)
  Lambda2Save=matrix(0,nmc,n)
  
  
  #Initialize
  Beta = y
  Tau = rgamma(1,c,rate= d^2)
  Sigma2 = 0.95*stats::var(y)
  Lambda2 = rexp(n,rate = Tau)
  
  
  for(t in 1:(nmc+burn)){
    
    #Update Beta
    f=Lambda2/(Lambda2+1)
    Beta=f*y+ sqrt(f*Sigma2)*rnorm(n)
    
    #Update Lambda
    Lambda2=1/rinvgauss(n,mean =sqrt(2*Sigma2*Tau/Beta^2),shape= 2*Tau)
    #Lambda2=1/rinv.gaussian(n, sqrt(2*Sigma2*Tau/Beta^2), lambda = 2*Tau)
    
    #Update Tau
    Tau=rgamma(1,shape= c+n,rate= sum(Lambda2)+d^2)
    
    
    #Update Sigma2
    Sigma2=1/rgamma(1, shape= n, rate= 0.5*(sum((y-Beta)^2)+sum(Beta^2/Lambda2)))
    
    
    
    
    
    #Save results
    if(t > burn){
      BetaSave[t-burn, ] = Beta
      TauSave[t-burn] = Tau
      Sigma2Save[t-burn] = Sigma2
      Lambda2Save[t-burn,]=Lambda2
      
    }
  }
  BetaHat = colMeans(BetaSave)
  BetaMedian = apply(BetaSave, 2, stats::median)
  TauHat = mean(TauSave)
  Sigma2Hat = mean(Sigma2Save)
  
  result <- list("BetaHat" = BetaHat, "BetaMedian" = BetaMedian,
                 "Sigma2Hat" = Sigma2Hat,
                 "TauHat" = TauHat, "BetaSamples" = BetaSave,
                 "TauSamples" = TauSave, "Sigma2Samples" = Sigma2Save,
                 "Lambda2Samples"=Lambda2Save)
  return(result)
  
  
}