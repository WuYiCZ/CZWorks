# average likelihood estimator for change point
al <- function(X){
  r <- dnorm(X, mean = 0, sd=1)/dnorm(X, mean = 1, sd=1)
  R <- c(0)
  for (i in 3:n) {
    R[2] <- r[1]
    R[i] <- R[i-1]*r[i-1]
  }
n <- length(X)
  k <- seq(1,n)
  k.al <- c()
  k.al <- sum(k*R)/sum(R)
  tau.al <- sample(c(floor(k.al),floor(k.al)+1), size=1, prob = c(1-k.al+floor(k.al),k.al-floor(k.al)))/n
 # tau.al <- sample(c(floor(k.al),floor(k.al)+1), size=1, prob = c((1-k.al+floor(k.al)),(k.al-floor(k.al))))/n
  return(tau.al)
}