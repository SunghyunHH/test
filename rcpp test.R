library(geoR) # inverse chi square를 위한 library
library(Rcpp)
#install.packages("RcppEigen")
#install.packages("RcppArmadillo")
library(RcppEigen)
library(RcppArmadillo)
#include <RcppEigen.h>
#include <Rcpp.h>
y <- list(A = c(62, 60, 63, 59),
          B = c(63, 67, 71, 64, 65, 66),
          C = c(68, 66, 71, 67, 68, 68),
          D = c(56, 62, 60, 61, 63, 64, 63, 59))

n.iter <- 3000
m <- 5  # num. of chains

# A,B,C,D로 theta가 4개, sample 3000개, 5번의 multiple chain이라
# array형태 (4 x 3000 x 5)
theta <- array(NA, dim=c(4, n.iter, m))
mu <- array(NA, dim=c(n.iter, m))
sigma2 <- array(NA, dim=c(n.iter, m))
tau2 <- array(NA, dim=c(n.iter, m))

# 각 chain에 intial value 설정
# unlist() -> list to vetor
for(int h=0, i<m; i++){
  # theta 와 mu는 data에서 sampling해서 초기값 줬음
  theta[,1,h] <- sample(unlist(y), 4)
  mu[1,h] <- sample(unlist(y), 1)
  sigma2[1,h] <- runif(1,1,10)
  tau2[1,h] <- runif(1,1,10)
}

# 각 chain별로 gibbs sampling 진행
for(h in 1:m){
  for(i in 2:n.iter){
    for(j in 1:4){
      var_hat <- 1/(1/tau2[i-1,h]+length(unlist(y[j]))/sigma2[i-1,h])
      theta_hat <- var_hat*(1/tau2[i-1,h]*mu[i-1,h]+length(unlist(y[j]))/sigma2[i-1,h]*mean(unlist(y[j])))
      theta[j,i,h] <- rnorm(1, theta_hat, sd=sqrt(var_hat))  
    }
    
    mu[i,h] <- rnorm(1, mean(theta[,i,h]), sd=sqrt(tau2[i-1,h]/4))
    
    sigma2_hat <- mean(unlist(sapply(1:4, function(x) unlist(y[x])-theta[x,i,h]))^2)
    sigma2[i,h] <- rinvchisq(1, length(unlist(y)), sigma2_hat)
    
    tau2_hat <- sum((theta[,i,h]-mu[i,h])^2)/(4-1)
    tau2[i,h] <- rinvchisq(1, 4-1, tau2_hat)
  }
}

plot(sqrt(sigma2)[,1], type="l", xlab="Iterations", ylab=expression(sigma^2))
for(t in 2:m){
  points(sqrt(sigma2)[,t], type="l", col=t)
  Sys.sleep(1)
}

# Theta_1
plot(theta[1, ,1], type="l", main=expression(theta[1]), ylab=expression(theta[1]), xlab="Iterations")
for(t in 2:m){
  points(theta[1,,t], type="l", col=t)
  Sys.sleep(1)
}

# Theta_4
plot(theta[4, ,1], type="l", main=expression(theta[4]), ylab=expression(theta[4]), xlab="Iterations")
for(t in 2:m){
  points(theta[4,,t], type="l", col=t)
  Sys.sleep(1)
}
