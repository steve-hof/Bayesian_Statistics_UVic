rm(list=ls())

library(coda)
library(R.matlab)
library(spate)
library(MASS)
library(statmod)
setwd("/Users/stevehof/Documents/School/Fall_2020/STAT_460/Assignments/Quality_Work")

update.beta = function(tau.curr, bet) {
  D.inv = solve((1 / tau.curr) * diag(p))
  A = partial.A + D.inv
  A.inv = solve(A)
  mu = A.inv %*% t(X) %*% y.tilde
  Sig = sig2 * A.inv
  mvrnorm(n = 1, mu = mu, Sigma = Sig)
}

update.tau = function(bet, lam, tau) {
  for(j in 1:p) {
    mu.prime = sqrt(lam^2 * sig2 / bet[j]^2)
    tau[j] = (rinvgauss(n = 1, mean = mu.prime, shape = lam^2))^2
  }
  return(1 / tau)
}

gibbs = function(n_iter, init, priors) {
  beta.out = matrix(data = NA, nrow = n_iter, ncol = p)
  beta.curr = init$beta
  beta.out[1, ] = beta.curr
  tau.curr = init$tau
  deb = 12
  
  for(k in 2:n_iter) {
    tau.curr = update.tau(bet = beta.curr, lam = priors$lam, tau = tau.curr)
    beta.curr = update.beta(tau.curr = tau.curr, bet = beta.curr)
    beta.out[k, ] = beta.curr
  }
  colnames(beta.out) = beta.names
  return(beta.out)
}

dat = readMat("Dataset_FinalProject_2.mat")
X = dat$X
y = dat$Y
sig2 = 1
n = length(y) # maybe move this and below to function?
p = dim(X)[2]
y.tilde = y - matrix(rep(1, n), nrow = n, ncol = 1)
partial.A = t(X) %*% X
beta.names = c("beta1", "beta2", "beta3", "beta4", "beta5",
               "beta6", "beta7", "beta8", "beta9", "beta10")

priors = list()
init = list()
n_iter = 10000
model = lm(y ~ X - 1)
init$beta = model$coefficients
init$tau = rep(1, p)

priors$sig2 = 1
priors$lam = p / sum(abs(init$beta))

post = gibbs(n_iter, init, priors)
saveRDS(post, file = "post_dec11.Rds")

calc.rhat = function(m, nchain, J, chain) {
  rhat = numeric(J)
  for (j in 1:J) {
    psi.mean = mean(chain[, j])
    psi.bar = numeric(m)
    aux.w = numeric(m)
    for (k in 1:m) {
      sub.chain = chain[seq((k - 1) * nchain + 1, k * nchain, 1), j]
      psi.bar[k] = mean(sub.chain)
      aux.w[k] = (1 / (nchain - 1)) * sum((sub.chain - mean(sub.chain))^2)
    }
    B = (nchain / (m - 1)) * (sum((psi.bar - psi.mean)^2))
    W = (1 / m) * sum(aux.w)
    VP = ((nchain - 1) / nchain) * W + (1 / nchain) * B
    rhat[j] = sqrt(VP / W)
    names(rhat) = beta.names
  }
  return(rhat)
}

rhat.post = calc.rhat(m = 5, nchain = 100, J = 10, chain = post)
rhat.post

model$coefficients  
post[10000, ]
reso = colMeans(post)
reso
coda::autocorr.plot(as.mcmc(post), lag.max = 30)
plot(as.mcmc(post))

