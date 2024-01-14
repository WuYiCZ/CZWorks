# RAK estimator for change point
RAK <- function(X){
  # X = x.sim[,1,1]
  Q <- matrix(NA,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (X[i] > X[j]) { Q[i,j] <- 1
      } else if (X[i] == X[j]) {
        Q[i,j] <- 0
      } else {
        Q[i,j] <- -1
      }
    }
  }
  RAK = rep(0,(n-1))
  for (k in 2:n) {
    RAK[k-1] <- abs(sum(Q[1:(k-1),k:n])) 
  }
  k.rak <- which.max(RAK)
  tau.rak <- k.rak/n
  return(tau.rak)
}
