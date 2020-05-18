library(truncdist)
library(extraDistr)
Kap.acceptanceprob=function(k.new,k.old,t,sig2,b,g){
  if(g==1){
    s2.new=(t^2+sig2*(1-2*k.new))/(2*k.new)
    s2.old=(t^2+sig2*(1-2*k.old))/(2*k.old)
    pi.new=log(dnorm(b,0,sd=sqrt(s2.new))*dbeta(k.new, 1/2,1))
    pi.old=log(dnorm(b,0,sd=sqrt(s2.old))*dbeta(k.old, 1/2,1))
  }else if(g==0){
    pi.new=log(dbeta(k.new, 1/2,1))
    pi.old=log(dbeta(k.old, 1/2,1))
  }
  
  q.old.new=log(dbeta(k.old,k.new,1))
  q.new.old=log(dbeta(k.new,k.old,1))
  ap=exp(min(0,pi.new+q.old.new-pi.old-q.new.old))
  
  return(ap)
  
}


sig.acceptanceprob=function(sig2.new,sig2.old,y.s,b,t,g,k){
  sum.new=-log(sig2.new)+dhcauchy(t-sqrt(sig2.new), sigma = sqrt(sig2.new), log = TRUE)
  sum.old=-log(sig2.old)+dhcauchy(t-sqrt(sig2.old), sigma = sqrt(sig2.old), log = TRUE)
  for (i in 1:length(y.s)){
    s2.new=(t^2+sig2.new*(1-2*k[i]))/(2*k[i])
    s2.old=(t^2+sig2.old*(1-2*k[i]))/(2*k[i])
    if(g[i]==1){
    sum.new=sum.new+log(dnorm(y.s[i],b[i],sqrt(sig2.new)))+
      log(dnorm(b[i],0,sqrt(s2.new)))
    sum.old=sum.old+log(dnorm(y.s[i],b[i],sqrt(sig2.old)))+
      log(dnorm(b[i],0,sqrt(s2.old)))
    }else if(g[i]==0){
      sum.new=sum.new+log(dnorm(y.s[i],0,sqrt(sig2.new)))
      sum.old=sum.old+log(dnorm(y.s[i],0,sqrt(sig2.old)))
    }
    }
  log.q.old.new=log(dtrunc(sig2.old, spec = "gamma", a=0, b=t^2, shape=sig2.new, 
                           rate=1))
  log.q.new.old=log(dtrunc(sig2.new, spec = "gamma", a=0, b=t^2, shape=sig2.old, 
                           rate=1))
  ap=exp(min(0,sum.new+log.q.old.new-sum.old-log.q.new.old))
  return (ap)
}


tau.acceptanceprob=function(t.new,t.old,k,sig2,b,g){
  log.pi.new=dhcauchy(t.new-sqrt(sig2),sigma = sqrt(sig2),log = TRUE)
  log.pi.old=dhcauchy(t.old-sqrt(sig2),sigma = sqrt(sig2),log = TRUE)
  for (i in 1:length(k)){
    if(g[i]==1){
    s2.new=(t.new^2+sig2*(1-2*k[i]))/(2*k[i])
    s2.old=(t.old^2+sig2*(1-2*k[i]))/(2*k[i])
    log.pi.new=log.pi.new+log(dnorm(b[i],0,sd=sqrt(s2.new)))
    log.pi.old=log.pi.old+log(dnorm(b[i],0,sd=sqrt(s2.old)))  
    }else if(g[i]==0){
      log.pi.new=log.pi.new+0
      log.pi.old=log.pi.old+0
    }
    
    }
  
  log.q.old.new= dhcauchy(t.old-sqrt(sig2),sigma = sqrt(sig2),log = TRUE)
  log.q.new.old= dhcauchy(t.new-sqrt(sig2),sigma = sqrt(sig2),log = TRUE)
  ap=exp(min(0,log.pi.new+log.q.old.new-log.pi.old-log.q.new.old))
  return (ap)
}



SB.normal.means <- function(y,burn=1000,nmc=5000,tau=3){
  
  n <- length(y)
  BetaSave = matrix(0, nmc, n)
  TauSave = rep(0, nmc)
  Sigma2Save = rep(0, nmc)
  GamSave=matrix(0,nmc,n)
  pSave = rep(0, nmc)
  
  #Initialize
  Beta = y
  Tau = tau
  Sigma2 = 0.5#0.95*stats::var(y)
  Lambda = 1/abs(y)^2
  Kappa=1/(1+Lambda)
  p=0.5
  Gam=rbern(n,prob = p)
  
  for(t in 1:(nmc+burn)){
    
    if (t%%1000 == 0){print(t)}
    
    #update block beta
    a=2*Kappa/(Tau^2+(1-2*Kappa)*Sigma2)
    b=Gam^2/Sigma2
    s=sqrt(1/(a+b))
    m=s^2*Gam*y/Sigma2
    for (ind in 1:n){
      if(Gam[ind]==0){
        Beta[ind]=0
      }else if (Gam[ind]==1){
        Beta[ind]=rnorm(1,mean = m[ind],sd=s[ind])
       }
    }
    
    #update kappa
    #Metropolis hastings with proposal Beta(kappa.prev,1) so that mean of proposal 
    #is (kappa.prev)/(1+kappa.prev)
    Kappa.old=Kappa
    for (ind in 1:n){
        Kap.new=rbeta(1,Kappa.old[ind],1)
        ap=Kap.acceptanceprob(Kap.new, Kappa.old[ind], Tau, Sigma2,Beta[ind],Gam[ind])
        u=runif(1)
        if (u<=ap){
          Kappa[ind]=Kap.new
        }else {
          Kappa[ind]=Kappa.old[ind]
        }
    }
    
    #function(g.new,g.old,dat,b,sig2,prob,k,t)
    #update Gam
    Gam.old=Gam
    
    for (ind in 1:n){
      stddev=sqrt((Tau^2+Sigma2)/(2*Kappa[ind]))
      p1=dnorm(y[ind],0,stddev)*p
      p0=(1-p)*dnorm(y[ind],0,sqrt(Sigma2))
      Gam.prob.1=p1/(p1+p0)
      Gam[ind]=rbern(1,prob = Gam.prob.1)
    }
    
    #Update p
    p=rbeta(1,1+sum(Gam),n+1-sum(Gam))
    
    
    #update sigma
    Sigma2.old=Sigma2
    Sigma2.new=rtrunc(1, spec = "gamma", a=0, b=Tau^2, shape=Sigma2.old, rate=1)#rgamma(1, Sigma2.old, rate = 1)
    #Sigma2.new=Sigma.new^2
    ap=sig.acceptanceprob(Sigma2.new,Sigma2.old, y, Beta, Tau, Gam,Kappa)
    u=runif(1)
    if (u<=ap){
      Sigma2=Sigma2.new
    }else {
      Sigma2=Sigma2.old
    }
    
    #update tau
    tau.old=Tau
    tau.new=rhcauchy(1, sqrt(Sigma2))+sqrt(Sigma2)
    ap=tau.acceptanceprob(tau.new,tau.old, Kappa, Sigma2,Beta,Gam)
    u=runif(1)
    if (u<=ap){
      Tau=tau.new
    }else {
      Tau=tau.old
    }
    
    #Save results
    if(t > burn){
      BetaSave[t-burn, ] = Beta
      TauSave[t-burn] = Tau
      Sigma2Save[t-burn] = Sigma2
      GamSave[t-burn, ]=Gam
      pSave[t-burn]=p
    }
    
  }#end MCMC
  
  #Compute results
  BetaHat = colMeans(BetaSave)
  BetaMedian = apply(BetaSave, 2, stats::median)
  TauHat = mean(TauSave)
  Sigma2Hat = mean(Sigma2Save)
  pip= apply(GamSave, 2, mean)
  
  result <- list("BetaHat" = BetaHat, "BetaMedian" = BetaMedian,
                 "Sigma2Hat" = Sigma2Hat,
                 "TauHat" = TauHat, "BetaSamples" = BetaSave,
                 "TauSamples" = TauSave, "Sigma2Samples" = Sigma2Save, 
                 "GamSamples"=GamSave,"PIP"=pip, "pSample"=pSave )
  return(result)
}