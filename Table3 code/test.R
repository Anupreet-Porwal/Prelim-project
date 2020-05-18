library(EbayesThresh)
library(extraDistr)
library(LaplacesDemon)
library(bayeslm)
library(mltools)
library(statmod)
library(truncdist)
library(data.table)

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

##### Data generation

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


data.generation <- function(p, xi,w,tau=3){
  gam=rbern(p,prob = w)
  theta=ifelse(gam, rst(p, mu=0, sigma = tau, nu=xi), 0)
  y=rnorm(p, theta, sd=1)
  data <- list("y"=y,"theta"=theta)
  return (data)
} 



#args = commandArgs(TRUE)
p <- 250
xi.c <- c(2,10)
w.c <- c(0.05,0.2,0.5)

n.sim <- 5

se_mat <- array(NA, dim = c(9,length(xi.c),length(w.c)))
methods <- c("MLE","Double-Exponential","NEG(c=4,d=3)","NEG(c=2,d=3)",
"NEG(c=1,d=3)","NEG(c=0.5,d=3)","Empirical Bayes",
"NEG(best fixed c,d)","Horseshoe")
dimnames(se_mat) <-  list(methods,c(paste("eps=",xi.c[1]),
                                     paste("eps=",xi.c[2])),
                           c(paste("w=",w.c[1]),paste("w=",w.c[2]),paste("w=",w.c[3])))
packages.list <- c('bayeslm','EbayesThresh', 'statmod', 'truncdist',
                   'extraDistr','LaplacesDemon','mltools')

cl<-detectCores()
registerDoParallel(cl)

for (i in 1:length(w.c)){
  for (j in 1:length(xi.c)){
    se_mat2 <- matrix(NA, nrow = n.sim, ncol = length(methods))
    colnames(se_mat2) <- methods
    
    results= foreach(sim =1:n.sim, .combine = cbind,
                     .packages = packages.list) %dopar% {
      
      if (sim%%1 == 0){print(paste("sim:",sim,", w=",w.c[i],",zeta=",xi.c[j]))}
      System$getHostname()
      # dat <- data.generation(p, xi=xi.c[j],w=w.c[i])
      # 
      # mle <- dat$y # Since X=I, theta_mle=y
      # 
      # de.mod <- quiet(bayeslm(dat$y, X=diag(p), prior = "laplace", icept = FALSE))
      # # de.mod2 <-  blasso(X=diag(p), y=dat$y, case="default", icept=FALSE)
      # 
      # hs.mod <- quiet(HSnormalmeans(dat$y))
      # 
      # eb.mod <- ebayesthresh(dat$y,verbose = TRUE)
      # 
      # neg.mod1 <- neg.normal.means(dat$y, c=4,d=3)
      # 
      # neg.mod2 <- neg.normal.means(dat$y, c=2,d=3)
      # 
      # neg.mod3 <- neg.normal.means(dat$y, c=1,d=3)
      # 
      # neg.mod4 <- neg.normal.means(dat$y, c=0.5,d=3)
      # 
      # 
      # c.seq=seq(0.5,8,0.5)
      # d.seq=seq(0.1,10,0.1)
      # se.curr=10000
      # for (c.s in 1:length(c.seq)){
      #   for(d.s in 1:length(d.seq)){
      #     neg.mod <- neg.normal.means(dat$y, c=c.seq[c.s],d=d.seq[d.s])
      #     se.new = mse(neg.mod$BetaHat,dat$theta)*p
      #     if(se.new<=se.curr){
      #       se.curr=se.new
      #     }
      #   }
      # }
      # 
      # # squared error loss 
      # temp_se=rep(0, 9)
      # temp_se[1] <-  mse(mle,dat$theta)*p
      # temp_se[2] <-  mse(de.mod$fitted.value,dat$theta)*p
      # temp_se[3] <-  mse(neg.mod1$BetaHat,dat$theta)*p
      # temp_se[4] <-  mse(neg.mod2$BetaHat,dat$theta)*p
      # temp_se[5] <-  mse(neg.mod3$BetaHat,dat$theta)*p
      # temp_se[6] <-  mse(neg.mod4$BetaHat,dat$theta)*p
      # temp_se[7] <-  mse(eb.mod$muhat,dat$theta)*p
      # temp_se[8] <-  se.curr
      # temp_se[9] <-  mse(hs.mod$BetaHat, dat$theta)*p
      # 
      # data.frame(sim=temp_se)
      # 
      
      # se_mat2[sim,1] <-  mse(mle,dat$theta)*p
      # se_mat2[sim,2] <-  mse(de.mod$fitted.value,dat$theta)*p
      # se_mat2[sim,3] <-  mse(neg.mod1$BetaHat,dat$theta)*p
      # se_mat2[sim,4] <-  mse(neg.mod2$BetaHat,dat$theta)*p
      # se_mat2[sim,5] <-  mse(neg.mod3$BetaHat,dat$theta)*p
      # se_mat2[sim,6] <-  mse(neg.mod4$BetaHat,dat$theta)*p
      # se_mat2[sim,7] <-  mse(eb.mod$muhat,dat$theta)*p
      # se_mat2[sim,8] <-  se.curr
      # se_mat2[sim,9] <-  mse(hs.mod$BetaHat, dat$theta)*p
    }
    # Insert values into table for given method and eps/w choice
    #se_mat[ ,j,i]=rowMeans(results) 
    # for (k in 1:length(methods)){
    #    se_mat[k,j,i]=mean(results[ ,k])
    #  }
  }
}

#fwrite(se_mat[ , ,1],file = "output1.csv")
#fwrite(se_mat[ , ,2],file = "output2.csv")
#fwrite(se_mat[ , ,3],file = "output3.csv")



stopCluster(cl)


