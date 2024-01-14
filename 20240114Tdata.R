## t distribution data
rm(list = ls()) # remove (almost) everything in the working environment.

### seting path
setwd("C:\\Users\\lenovo\\Desktop\\JoAS2\\Pro")
### loding libraries ...
library(MASS)
library(mixAK)

### loding change point eatimation algorithm 
source("ISA.R")
source("CUSUM.R")
source("RAK.R")
source("PCUS.R")

### generating data
N <- 1000; n <- 100; m <- 5
theta <- 0.6
v = 4 # degrees of freedom of the multivariate t distribution
# the ratio \tau of change point
tau1 <- c(0,0.25,0.5,0.75,1)
m <- length(tau1)
delta <- 0.1 ### amounts of structural change
Sigma <- matrix(0, n+m, n+m) ### set co-variances matrix 
for (i in 1:(m+n-1)) {
  Sigma[i,i+1] <- -theta
}
for (i in 2:(m+n)) {
  Sigma[i,i-1] <- -theta
}
for (i in 1:(m+n)) {
  Sigma[i,i] <- 1 + theta^2
}

mu1 <- 2; mu2 <- 2 + n^(-delta)

x.sim <- array(NA, dim = c(n,N,m))
for (t1 in 1:5) {
  k <- floor(n*tau1[t1]) # location of change point
  mu <- rep(0, (m+n))
  e <- rMVT(N,v,mu,1,Sigma)$x # N row m+n col
  epsilon <- rep(0, (n+m))
  x <- matrix(NA, n, N)
  for (j in 1:N) {
    for (i in 1:n) {
      epsilon[i] <- e[j, i:i+m-1]
    }
    varepsilon <- runif(n, 0, .5)
    Y <- rep(0,n)
    for (i in 1:n) {
      Y[i+1] = Y[i]*varepsilon[i] + epsilon[i]
    }
    Y <- Y[-1] # the m-ANA random variables
    x[1:k, j] <- Y[1:k] + mu1
    if (k < n) {
      x[(k + 1):n, j] <- Y[(k + 1):n] + mu2
    } 
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
print(table.wspm)

#################### CUSUM estimator #########################
Abias.cusum <- MSE.cusum <- c(NA,5)
alpha0 <- 0
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
Abias.RAK <- MSE.RAK <- c(NA, m)
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
### getting means
