# average likelihood estimator for change point
den.f <- function(h,t,e){
  res = 1/(sqrt(2*pi)*h)*mean(exp(-(t-e)^2/(2*h^2))) 
  return(res)
}

al <- function(X){
 # X <- x.sim[,2,3]
  n <- length(X)
  hn <- n^(-1/5);hm <- n^(-1/4)
  for (i in 1:n) {
    r1[i] <- den.f(h=hn,t=X[i],e=X[1:(n/2)]);
    r2[i] <- den.f(h=hm,t=X[i],e=X[(n/2+1):n]);
  }
  r <- r1/r2
  # r <- dnorm(X, mean = 0, sd=1)/dnorm(X, mean = 1, sd=1)
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
  tau.al
  # tau.al <- sample(c(floor(k.al),floor(k.al)+1), size=1, prob = c((1-k.al+floor(k.al)),(k.al-floor(k.al))))/n
  return(tau.al)
}