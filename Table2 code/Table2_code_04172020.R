library(horseshoe)

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
# 11 seed works the best with 041720 code with tau=1 with burn=3000; 5000 is still not working
#set.seed(11)
n <- c(25,50,100,200,500,1000,2000,5000)#,10000)
weights.table <- matrix(0,nrow = length(n),ncol = 10)
FP.table <- matrix(0,nrow = length(n),ncol = 1)
weights.table.SB <- matrix(0,nrow = length(n),ncol = 10)
FP.table.SB <- matrix(0,nrow = length(n),ncol = 1)


for (i in 1:length(n)){
  print(paste0("Noises:",n[i]))
  truth <- c(rep(0, n[i]), seq(0.5,5,0.5))
  data <- truth + rnorm((n[i]+10))
  
  # Note that this is using manually written MCMC sampler for Horseshoe 
  res.HS2 <- HSnormalmeans(data)#quiet(HS.normal.means(data, method.tau = "truncatedCauchy", method.sigma = "Jeffreys"))
  #Plot the posterior mean against the data (signals in blue)
  #plot(data, res.HS2$BetaHat, col = c(rep("black", n[i]), rep("blue", 10)))
  #Find the selected betas (ideally, the last 20 are equal to 1)
  var.sel <-  HS.var.select(res.HS2, data, method = "threshold",threshold = 0.5)
  proportion <- abs(res.HS2$BetaHat/data)
  weights.table[i, ] <- proportion[(n[i]+1):(n[i]+10)]
  FP.table[i,1] <- sum(proportion[1:n[i]]>=0.5)#sum(var.sel[1:n[i]]>0) 
  res.SB <- SB.normal.means(data,burn = 3000,tau=1)
  plot(res.SB$Sigma2Samples,type="l")
  print(length(unique(res.SB$Sigma2Samples)))
  weights.table.SB[i, ] <- res.SB$PIP[(n[i]+1):(n[i]+10)]
  FP.table.SB[i,1] <- sum(res.SB$PIP[1:n[i]]>=0.5)
  
 }

sweights.HS <-  round(weights.table*100)
pip.SB <- round(weights.table.SB*100)

colnames(sweights.HS)=colnames(pip.SB)=seq(0.5,5,0.5)
rownames(sweights.HS)=rownames(pip.SB)=n

rownames(FP.table)=rownames(FP.table.SB)=n
colnames(FP.table)=colnames(FP.table.SB)="FP"

summary.HS <- cbind(sweights.HS,FP.table)
summary.SB <- cbind(pip.SB,FP.table.SB)

# > xtable(summary.HS, digits = 0)
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Fri Apr 17 21:38:41 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrrr}
# \hline
# & 0.5 & 1 & 1.5 & 2 & 2.5 & 3 & 3.5 & 4 & 4.5 & 5 & FP \\ 
# \hline
# 25 & 15 & 12 & 18 & 22 & 15 & 16 & 24 & 32 & 11 & 41 & 0 \\ 
# 50 & 27 & 23 & 16 & 25 & 33 & 57 & 85 & 93 & 93 & 96 & 0 \\ 
# 100 & 9 & 8 & 9 & 26 & 24 & 30 & 68 & 41 & 91 & 18 & 0 \\ 
# 200 & 13 & 14 & 21 & 78 & 61 & 92 & 80 & 70 & 92 & 95 & 1 \\ 
# 500 & 7 & 3 & 2 & 5 & 5 & 8 & 81 & 13 & 94 & 91 & 0 \\ 
# 1000 & 3 & 2 & 3 & 2 & 5 & 13 & 42 & 76 & 85 & 88 & 2 \\ 
# 2000 & 2 & 2 & 1 & 46 & 3 & 18 & 8 & 13 & 93 & 92 & 1 \\ 
# 5000 & 0 & 2 & 0 & 1 & 2 & 4 & 5 & 6 & 89 & 90 & 6 \\ 
# \hline
# \end{tabular}
# \end{table}
# > xtable(summary.SB, digits = 0)
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Fri Apr 17 21:39:34 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrrr}
# \hline
# & 0.5 & 1 & 1.5 & 2 & 2.5 & 3 & 3.5 & 4 & 4.5 & 5 & FP \\
# \hline
# 25 & 7 & 6 & 8 & 14 & 8 & 9 & 15 & 24 & 5 & 33 & 0 \\
# 50 & 21 & 15 & 11 & 20 & 29 & 59 & 94 & 100 & 100 & 100 & 0 \\
# 100 & 5 & 4 & 5 & 19 & 16 & 25 & 80 & 37 & 98 & 10 & 0 \\
# 200 & 6 & 8 & 11 & 86 & 61 & 99 & 84 & 75 & 99 & 100 & 1 \\
# 500 & 6 & 2 & 1 & 3 & 2 & 6 & 85 & 10 & 100 & 99 & 0 \\
# 1000 & 2 & 1 & 1 & 1 & 2 & 8 & 39 & 67 & 86 & 99 & 0 \\
# 2000 & 2 & 1 & 1 & 22 & 1 & 9 & 4 & 12 & 100 & 100 & 1 \\
# 5000 & 14 & 95 & 31 & 76 & 73 & 99 & 98 & 100 & 100 & 100 & 440 \\
# \hline
# \end{tabular}
# \end{table}