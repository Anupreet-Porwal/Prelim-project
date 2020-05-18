library(kernlab)

generate_kernel_data <- function(t){
  return(sin(t)/t+rnorm(length(t),mean = 0,sd=0.15))
}

seed(1)
sim=100
sse.rvm=sse.HS=count= 0

for (i in 1:sim){
  t <- runif(100,min = -20,max = 20)
  yt <- generate_kernel_data(t) 
  foo <- rvm(t,yt)
  
  sig <- kpar(kernelf(foo))$sigma
  rbf <- rbfdot(sigma=sig)
  xmat <- kernelMatrix(rbf,t)
  
  
  HS.mod <- try(log("a"),silent = TRUE)
  while(class(HS.mod)=="try-error"){
  HS.mod <- try(HS.regression(X=xmat,y=yt))
  }
  
  ttest <- runif(100,min = -20,max = 20)
  y.test <- sin(ttest)/ttest
  
  rvm.pred <- predict(foo,ttest)
  xmat.test <- kernelMatrix(rbf, ttest,t)
  HS.pred <- xmat.test %*% HS.mod$BetaHat
  
  sse.HS <- sse.HS+sum((y.test-HS.pred)^2)
  sse.rvm <- sse.rvm+sum((y.test-rvm.pred)^2)
  if(sum((y.test-HS.pred)^2)<=sum((y.test-rvm.pred)^2)){count =count +1}
  
  
  yt.rvm <- predict(foo, sort(t))
  yt.hs <- xmat %*% HS.mod$BetaHat
  plot(t, yt, type ="p",xlab="x",ylab="sin(x)/x",pch=20)
  lines(sort(t), yt.rvm, col="red",type = "l")
  lines(sort(t), yt.hs[order(t)], col="blue",type = "l")
  lines(sort(t),sin(sort(t))/sort(t))
  legend("topleft", legend=c("Truth", "Horseshoe (SSE=0.91)","RVM (SSE=3.24)"),
         col=c("black","blue", "red"), lty=1, cex=0.8)
  # 
  
}



# > count
# [1] 58
# > sse.HS
# [1] 77.95661
# > sse.rvm
# [1] 93.27915