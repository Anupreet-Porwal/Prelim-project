library(huge)
library(CVglasso)
library(gRim)
library(mixggm)
library(sparsebn)

# Need to run this again to get output 

setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Spring 2020/STAT572/code")
# L <- huge.generator(n = 86, d = 10, graph = "hub", g = 4)
# data <- L$data
# L <- data("cytometryContinuous")
# dat <- sparsebnData(cytometryContinuous$data, type = "c", ivn = cytometryContinuous$ivn)
# data <- dat$data
data <- read.table("Vanguard data.txt")


train <- as.matrix(data[1:60,])
p <- ncol(train)
test <- data[61:86,]

# change made here!! uncomment to reach back to normal
pr.mat <- solve(t(train)%*%train)#/nrow(train))
mle.prec.est <- 0.5*(pr.mat+t(pr.mat))



#write.table(train,"train.txt",col.names = FALSE,row.names = FALSE)

lasso.and <- huge(train, method = "mb",sym="and")
out.and <- huge.select(lasso.and)
lasso.and.f <- huge.mb(train, lambda = out.and$opt.lambda,sym="and")
and.path.est <- as(lasso.and.f$path[[1]],"matrix")
and.est <- as(lasso.and.f$beta[[1]],"matrix")#*as(lasso.and.f$path[[1]],"matrix")
and.ggm.est <- fitGGM(S=solve(mle.prec.est),N=nrow(train), graph = and.path.est, model = "concentration")
and.prec.est <- and.ggm.est$omega


lasso.or <- huge(train, method = "mb",sym="or")
out.or <- huge.select(lasso.or)
lasso.or.f <- huge.mb(train, lambda = out.or$opt.lambda,sym="or")
or.est <- as(lasso.or.f$beta[[1]],"matrix")#*as(lasso.or.f$path[[1]],"matrix")
or.path.est <- as(lasso.or.f$path[[1]],"matrix")
or.ggm.est <- fitGGM(S=solve(mle.prec.est),N=nrow(train), graph = or.path.est, model = "concentration")
or.prec.est <- or.ggm.est$omega


glasso.mod <- CVglasso(X=train)
glasso.prec.est <- glasso.mod$Omega


T=diag(p)
D=rep(0, p)
D[1]=var(train[ ,1])
for (i in 2:p){
  y=train[ , i]
  x=as.matrix(train[ , 1:(i-1)])
  HS.mod <- HS.regression(x,y)
  T[i, 1:(i-1)]=-HS.mod$BetaHat
  D[i]=HS.mod$Sigma2Hat
}
HS.prec.est <- t(T)%*%diag(1/D)%*%T




test.sets <- combn(1:p,3)

se.sum <- function(test,t,a,lt,meth){
  y.act <- test[t,]
  y.pred <- y.act
  y.pred[a] <- NA
  
  gam.aa=lt[a,a]
  gam.ab=lt[a,-a]
  y.pred[a]=-solve(gam.aa)%*%gam.ab%*%t(as.matrix(y.pred[-a]))
  se=sum((y.act[a]-y.pred[a])^2)
    
  return (se)
}

ae.sum <- function(test,t,a,lt,meth){
  y.act <- test[t,]
  y.pred <- y.act
  y.pred[a] <- NA
  
  gam.aa=lt[a,a]
  gam.ab=lt[a,-a]
  y.pred[a]=-solve(gam.aa)%*%gam.ab%*%t(as.matrix(y.pred[-a]))
  ae=sum(abs(y.act[a]-y.pred[a]))
  
  return (ae)
}



lt.list <- list("mle"= mle.prec.est,"or"=or.prec.est,"and"=and.prec.est,
                "glasso"=glasso.prec.est, "HS"=HS.prec.est)#,"or"=or.est,"and"=and.est,"glasso"=glass.unichol,"HS"=HS.unichol)

#pred <- array(NA,dim=c(nrow(test.sets),ncol(test.sets),ncol(test)))
se.mat=rep(0,length(lt.list))
ae.mat=rep(0,length(lt.list))

for (k in 1:nrow(test)){
  print(k)
  for( j in 1:ncol(test.sets)){
    for (q in 1:length(se.mat)){
      se.mat[q]=se.mat[q]+se.sum(test,k,test.sets[ ,j],lt.list[[q]],meth = names(lt.list)[q])
      ae.mat[q]=ae.mat[q]+ae.sum(test,k,test.sets[ ,j],lt.list[[q]],meth = names(lt.list)[q])
    }
  }
}

output.se <- as.numeric(se.mat)/nrow(test.sets)
output.ae <- as.numeric(ae.mat)/nrow(test.sets)

# Order is MLE, OR, AND, Glasso, Horseshoe 
# > output.ae
# [1] 17599.562  6719.430  6379.122  6393.022  5537.986
# > output.se
# [1] 1376.6313  214.9977  182.4157  170.8892  155.7959
# > output.ae/output.ae[5]
# [1] 3.177972 1.213335 1.151885 1.154395 1.000000
# > output.se/output.se[5]
# [1] 8.836122 1.379996 1.170864 1.096879 1.000000
