# CUSUM estimator for change point
CUSUM <- function(X, alpha){
  # X=x.sim[,1,1]; alpha = 1/2
  n <- length(X) # sample size
k1 <- 1; k2 <- n - 1
  U <- rep(0, k2-k1+1 )
  for (k in  k1:k2) {
    U[k] <- (k*(n-k)/n)^(1-alpha)*(mean(X[1:k]) - mean(X[(k+1):n]))
  }
  U <- abs(U)
  kbo <- which.max(U); taubo <- kbo/n
  return(taubo)
}
