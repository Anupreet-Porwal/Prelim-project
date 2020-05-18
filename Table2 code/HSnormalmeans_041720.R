library(truncdist)

# This paper implements the Horeshoe sampler discussed by 
# Makalic, Enes, and Daniel F. Schmidt. "A simple sampler for the horseshoe estimator." 
# IEEE Signal Processing Letters 23.1 (2015): 179-182.
# It uses the relation between Cauchy and Inverse gamma to make conditional distributions
# of all the parameters involved conjugate 
# Another approach could be to do slice sampling which is implemented in another file
# HSnormalmean_041720_vAlt

HSnormalmeans <- function(y,burn=1000,nmc=5000,tau=1){
  
  n <- length(y)
  BetaSave = matrix(0, nmc, n)
  LambdaSave = matrix(0, nmc, n)
  TauSave = rep(0, nmc)
  Sigma2Save = rep(0, nmc)
  
  
  
  #Initialize
  Beta = y
  Tau = tau
  Sigma2 = 0.95*stats::var(y)
  Lambda = 1/abs(y)^2
  nu=1/rgamma(n, shape=1/2,rate=1)
  chi=1/rgamma(1,shape=1/2,rate=1)

  for(t in 1:(nmc+burn)){
    
    if (t%%1000 == 0){print(t)}
    
    # Update Beta
    num=(Lambda^2)*(Tau^2)
    f=num/(1+num)
    Beta=rnorm(n, mean = f*y, sd=sqrt(f*Sigma2))
    
    
    #Update Sigma2
    
    Sigma2.inv=rgamma(1,n,rate = 0.5*(sum((y-Beta)^2)+sum((Beta/Lambda)^2)/Tau^2))
    Sigma2=1/Sigma2.inv
    
    # Update Lambda
    b1=1/nu+Beta^2/(2*Tau^2*Sigma2)
    Lambda=sqrt(1/rgamma(n, shape = 1, rate = b1))
    
    # Update Tau
    Theta=Beta/Lambda
    b2=1/chi+sum(Theta^2)/(2*Sigma2)
    Tau=sqrt(1/rgamma(1,(n+1)/2,b2))
    
    # Update nu 
    nu=1/rgamma(n,1,1+1/Lambda^2)
    
    # Update chi
    chi=1/rgamma(1,1,1+1/Tau^2)
    
    
    #Save results
    if(t > burn){
      BetaSave[t-burn, ] = Beta
      TauSave[t-burn] = Tau
      Sigma2Save[t-burn] = Sigma2
      LambdaSave[t-burn, ] = Lambda
      
    }
  }
  BetaHat = colMeans(BetaSave)
  BetaMedian = apply(BetaSave, 2, stats::median)
  TauHat = mean(TauSave)
  Sigma2Hat = mean(Sigma2Save)
  
  result <- list("BetaHat" = BetaHat, "BetaMedian" = BetaMedian,
                 "Sigma2Hat" = Sigma2Hat,
                 "TauHat" = TauHat, "BetaSamples" = BetaSave,
                 "TauSamples" = TauSave, "Sigma2Samples" = Sigma2Save, "LambdaSamples"=LambdaSave)
  return(result)
  
}
