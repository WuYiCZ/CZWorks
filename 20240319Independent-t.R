## normaldata
rm(list = ls()) # remove (almost) everything in the working environment.

### seting path
setwd("C:\\Users\\lenovo\\Desktop\\二修\\Pro")

### loding change point estimation algorithm 
source("ISA.R")
source("CUSUM.R")
source("RAK.R")
source("PCUS.R")
source("ALt.R")

### generating data
N <- 1000; n <- 500; tau <- 0.25
tau1 <- c(0,0.25,0.5,0.75,1)
m <- length(tau1)
x.sim <- array(NA, dim = c(n,N,m))
for (t1 in 1:5) {
  k <- floor(n*tau1[t1]) # location of change point
  x <- matrix(NA, n, N)
  for (j in 1:N) {
    x1 <- rt(k, df=1);
    x2 <- rt(n-k,  df=1 + n^(-0.1)) 
    x[,j] <- c(x1,x2)
  }
  x.sim[,,t1] <- x
}

#################### WSPM estimator #########################
Abias.wspm <- MSE.wspm <- c(NA,5)
r0 <- 0.1 ### WSPM 
for (t1 in 1:5) {
  tauNhat <- c()
  for (t in 1:N) {
    tauNhat[t] <- ISA(X = x.sim[, t, t1], eta = 1, r=r0) 
  }
  # tauNha <- ISA(X = x.sim[, 1,2], eta = 1, r=2)
  # tauNha
  ### getting Abias and MSE for estimator 
  Abias.wspm[t1] <- mean(abs(tauNhat - tau1[t1]))
  MSE.wspm[t1] <- mean((tauNhat - tau1[t1])^2)
  print(t1)
}

table.wspm <- c()
for (t1 in 1:5) {
  table.wspm[2*t1 - 1] <- Abias.wspm[t1]
  table.wspm[2*t1] <- MSE.wspm[t1]
}
table.wspm <- round(table.wspm,4)
print(table.wspm)  # the 1'th row for r=0.1 in table 1;the 2'th row for r=1.5 in table 1


#################### average likelihood estimator #########################
Abias.al <- MSE.al <- c(NA,5)
for (t1 in 1:5) {
  tauNhat <- c()
  for (t in 1:N) {
    tauNhat[t] <- al(X = x.sim[, t, t1]) 
  }
  ### getting Abias and MSE for estimator 
  Abias.al[t1] <- mean(abs(tauNhat - tau1[t1]))
  MSE.al[t1] <- mean((tauNhat - tau1[t1])^2)
  print(t1)
}
t1
#plot(seq(1:200),c(x1,x2))

table.al <- c()
for (t1 in 1:5) {
  table.al[2*t1 - 1] <- Abias.al[t1]
  table.al[2*t1] <- MSE.al[t1]
}
table.al <- round(table.al, 4)
print(table.al)


#################### CUSUM estimator #########################
Abias.cusum <- MSE.cusum <- c(NA,5)
alpha0 <- 0.5
for (t1 in 1:5) {
  tauNhat <- c()
  # t1 = 2
  for (t in 1:N) {
    tauNhat[t] <- CUSUM(X = x.sim[, t, t1], alpha = alpha0) 
  }
  ### getting Abias and MSE for cusum estimator 
  Abias.cusum[t1] <- mean(abs(tauNhat - tau1[t1]))
  MSE.cusum[t1] <- mean((tauNhat - tau1[t1])^2)
  print(t1)
}

table.cusum <- c()
for (t1 in 1:5) {
  table.cusum[2*t1 - 1] <- Abias.cusum[t1]
  table.cusum[2*t1] <- MSE.cusum[t1]
}
table.cusum <- round(table.cusum, 4)
print(table.cusum)

#################### RAK estimator #########################
Abias.RAK <- MSE.RAK <- c(NA,5)
for (t1 in 1:5) {
  tauNhat <- c()
  # t1 = 2
  for (t in 1:N) {
    tauNhat[t] <- RAK(X = x.sim[, t, t1]) 
  }
  ### getting Abias and MSE for cusum estimator 
  Abias.RAK[t1] <- mean(abs(tauNhat - tau1[t1]))
  MSE.RAK[t1] <- mean((tauNhat - tau1[t1])^2)
  print(t1)
}

table.RAK <- c()
for (t1 in 1:5) {
  table.RAK[2*t1 - 1] <- Abias.RAK[t1]
  table.RAK[2*t1] <- MSE.RAK[t1]
}
table.RAK <- round(table.RAK, 4)
print(table.RAK)

#################### PCUS estimator #########################
Abias.PCUS <- MSE.PCUS <- c(NA,5)
for (t1 in 1:5) {
  tauNhat <- c()
  # t1 = 2
  for (t in 1:N) {
    tauNhat[t] <- PCUS(X = x.sim[, t, t1]) 
  }
  ### getting Abias and MSE for PCUS estimator 
  Abias.PCUS[t1] <- mean(abs(tauNhat - tau1[t1]))
  MSE.PCUS[t1] <- mean((tauNhat - tau1[t1])^2)
  print(t1)
}

table.PCUS <- c()
for (t1 in 1:5) {
  table.PCUS[2*t1 - 1] <- Abias.PCUS[t1]
  table.PCUS[2*t1] <- MSE.PCUS[t1]
}
table.PCUS <- round(table.PCUS, 4)
print(table.PCUS)
