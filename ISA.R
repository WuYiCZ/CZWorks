# Iterative searching algorithm

# See The weighted sum of powers in mean for eatimating a change point in linear processes
# with random coefficients for information about this R code. 

# Copyright (C) 2024 Shipeng Wu.

ISA = function (X, eta = 1, r) {
  if (r == 1 | r <= 0) { print("the r does not meet the requirements of the algorithm")}
  # X <- x.sim[,2,5]
  n <- length(X)
  nabla <- log(log(n))
  M <- max(abs(X)) + 1 ## step 1
  X <- X + M
  i <- 1; j <- 0; l <- n; dk <- 0
  muhat <- mean(X)
  khat <- c()
  if (r > 1) { # the case of r >1 
    while (dk < eta) {
      k1 <- max(floor(j + nabla),1); k2 <- min(ceiling(j + l - nabla),(n-1)) 
      w1 <- rep(-Inf, n)
      for (k in k1:k2) {
        mu1hat <- mean(X[1:k]); mu2hat <- mean(X[(k+1):n]); 
        w1[k] <-  k*mu1hat^r + (n-k)*mu2hat^r - n*muhat^r
      }
      which.max(w1)
      khat[i] <- which.max(w1)
      k3 <- max(floor(khat[i] - l^{1/2}),1)
      k4 <- min(ceiling(khat[i] + l^{1/2}), (n-1))
      lk34 <- length(k4-k3)
      w2 <- rep(-Inf, lk34)
      # k in (khat[i] - l^{1/2}, khat[i] + l^{1/2})) but k != khat[i]
      k34.range <- c(k3:k4)
      k34.range <- k34.range[!k34.range %in% c(khat[i])] 
      for (k in k34.range) {
        mu1hat <- mean(X[1:k]); mu2hat <- mean(X[(k+1):n]) 
        w2[k] <-  k*mu1hat^r + (n-k)*mu2hat^r - n*muhat^r
      }
      khat[i+1] <- which.max(w2)
      l <- 2*l^(1/2); j = khat[i] - l^(1/2)
      dk <- abs(khat[i+1] - khat[i])
      i <- i + 1
    }
    knhat <- khat[i-1]
  } else { # the case of 0 < r <1
    while (dk < eta) {
      k1 <- max(floor(j + nabla),1); k2 <- min(ceiling( j + l - nabla),(n-1)) 
      w1 <- rep(Inf,n)
      muhat <- mean(X) 
      for (k in k1:k2) {
        mu1hat <- mean(X[1:k]); mu2hat <- mean(X[(k+1):n]); 
        w1[k] <-  k*mu1hat^r + (n-k)*mu2hat^r - n*muhat^r
      }
      which.min(w1)
      khat[i] <- which.min(w1)
      k3 <- max(floor(khat[i] - l^{1/2}),1) 
      k4 <- min(ceiling(khat[i] + l^{1/2}),(n-1))
      lk34 <- length(k4-k3)
      w2 <- rep(Inf,n)
      k34.range <- c(k3:k4)
      k34.range <- k34.range[!k34.range %in% c(khat[i])]
      for (k in k34.range) {
        mu1hat <- mean(X[1:k]); mu2hat <- mean(X[(k+1):n]) 
        w2[k] <-  k*mu1hat^r + (n-k)*mu2hat^r - n*muhat^r
      }
      khat[i+1] <- which.min(w2)
      l <- 2*l^{1/2}; j = khat[i] - l^{1/2}
      dk <- abs(khat[i+1] - khat[i])
      i <- i + 1
    }
    knhat <- khat[i-1]
  }
  taunhat <- knhat/n
  return(taunhat)
}
# ISA(X=x.sim[,1,1],eta = 1, r=2)
