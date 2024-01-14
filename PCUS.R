# PCUS estimator for change point
PCUS <- function(X){
  # X = x.sim[,1,1] 
  n <- length(X)
  xbar <- mean(X)
  pcus <- rep(0,(n-1))
  for (k in 2:n) {
    pcus[k-1] <- ((k-1)*xbar - sum(X[1:(k-1)]))/sqrt((k-1)*(1-(k-1)/n))
  }
  k.pcus <- which.max(pcus) + 1
  tau.pcus <- k.pcus/n
  return(tau.pcus)
}
