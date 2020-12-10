rm(list=ls())

library(coda)
library(R.matlab)
library(spate)
library(MASS)
library(statmod)
setwd("/Users/stevehof/Documents/School/Fall_2020/STAT_460/Assignments/Quality_Work")


dat = readMat("Dataset_FinalProject_2.mat")

X = dat$X
y = dat$Y
n = length(y)
p = dim(X)[2]
sig2 = 1
tau = rep(1, p)
partial.A = t(X) %*% X

indx_gen = function(j, M) {
  n = dim(M)[2]
  if (j == 1) return(M[1, ])
  if (j == n) return(M[2, ])
  return(c(M[2, 1 : j-1], M[1, j:n]))
}

update.beta = function(tau.curr, bet) {
  D = (1 / tau.curr) * diag(p)
  D.inv = solve(D)
  A = partial.A + D.inv
  deb = 12
  y.tilde = numeric(n)
  for(i in 1:n) {
    for(j in 1:p) {
      y.tilde[i] = (y[i] - sum(setdiff(X[i,], X[i,j]) * setdiff(bet, bet[j])))  
      deb = 12
    }
    deb = 12
  }
  A.inv = solve(A)
  y.tilde.mat = matrix(NA, nrow = 300, ncol = 1)
  y.tilde.mat[, 1] = y.tilde
  mu = A.inv %*% t(X) %*% y.tilde.mat
  sd = sig2 * A.inv
  mvrnorm(n = 1, mu, sd)
}

update.tau = function(bet, lam, tau) {
  lam = 1
  for(j in 1:p) {
    mu.prime = sqrt(lam^2 * sig2 / bet[j]^2)
    tau[j] = rinvgauss(n = 1, mean = mu.prime, shape = lam)
  }
  return(1 / tau)
}


model = lm(y ~ X - 1)
beta.init = model$coefficients
tau = update.tau(beta.init, lam = 1, tau = tau)
bet = update.beta(tau.curr = tau, bet = beta.init)
deb = 12

  
