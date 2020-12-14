rm(list=ls())

library(coda)
library(R.matlab)
library(spate)
setwd("/Users/stevehof/Documents/School/Fall_2020/STAT_460/Assignments/")

dat = readMat("Quality_Work/Dataset_FinalProject_2.mat")
dat2 = readMat("not_include/mik_data.mat")
X = dat2$X
y = dat2$Y

update_pi = function(delta_j, a_pi, b_pi) {
  shape1 = sum(delta_j) + a_pi/2
  shape2 = sum(1-delta_j) + b_pi/2
  rbeta(n = 1, shape1, shape2)
}

update_delta = function(beta_j, pi, tau, eps) {
  p0 = ((1 - pi) * tau / sqrt(eps)) * exp(-1/(2*eps) * beta_j^2)
  p1 = pi * exp(-1/(2*tau^2) * beta_j^2)
  rob = p1 / (p0 + p1)
  bobbyd = rbinom(n = 1, size = 1, prob = rob)
  return(bobbyd)
}

update_beta = function(tau, eps, X, y, delta_j, j, bet, mutop = 0, mubot = 0) {
  partial.mutop = mutop
  partial.mubot = mubot
  for(i in 1:dim(X)[1]) {
    # partial.mutop = partial.mutop + X[i, j] * (y[i] - sum(setdiff(X[i,], X[i,j]) * 
    #                                                         setdiff(bet, bet[j])))
    partial.mutop = partial.mutop + X[i, j] * (y[i] - sum(setdiff(X[i,], X[i,j]) * 
                                                            setdiff(bet, bet[j])))
    partial.mubot = partial.mubot + X[i, j]^2
    deb = 12
  }
  deb = 12
  if(delta_j == 0) {
    sd = 1 + eps * partial.mubot
    mu = eps * partial.mutop / sd
    bob = (rnorm(n = 1, mean = mu, sd = sqrt(sd^-1)))
    return(bob)
  }
  if(delta_j == 1) {
    sd = 1 + tau^2 * partial.mubot
    mu = tau^2 * partial.mutop / sd
    bob = (rnorm(n = 1, mean = mu, sd = sqrt(sd^-1)))
    return(bob)
  }
}

gibbs = function(y, X, n_iter, init, priors) {
  # initialize
  bob = dim(X)[2]
  beta.out = matrix(data = NA, nrow = n_iter, ncol = dim(X)[2])

  delta.curr = init$delta_j
  beta.curr = init$beta
  beta.out[1, ] = beta.curr
  deb = 12
  # Gibbs sampler
  for(k in 2:n_iter) {
    pi.curr = update_pi(delta_j = delta.curr, a_pi = priors$a_pi, 
                        b_pi = priors$b_pi)
    
    for(j in 1:length(beta.curr)) {
      delta.curr = update_delta(beta_j = beta.curr[j], pi = pi.curr, 
                                tau = priors$tau, eps = priors$eps)
      # TODO: I think my issue is that I'm doing all the betas at once or something
      #       because my traceplot is pretty darn smooth which doesn't make a lot of sense.
      beta.curr[j] = update_beta(tau = priors$tau, eps = priors$eps, X = X, y = y,
                                delta_j = delta.curr, j = j, bet = beta.out[k-1])
      deb = 12
      beta.out[k, j] = beta.curr[j]
    }
    # beta.out[k, ] = beta.curr
  }
  return(beta.out)
}

# Work
priors = list()
init = list()
n_iter = 10000

model = lm(y ~ X - 1)
summary(model)
init$beta = model$coefficients
init$delta_j = 0

priors$a_pi = 1
priors$b_pi = 1
priors$tau = 10
priors$eps = 10^-4
set.seed(53)
post = gibbs(y = y, X = X, n_iter = n_iter, init = init, priors = priors)
colnames(post) = c("beta1", "beta2", "beta3", "beta4", "beta5", 
                    "beta6", "beta7", "beta8", "beta9", "beta10")

tail(post)
# plot(as.mcmc(post))
# summary(as.mcmc(post))
# coda::autocorr.plot(as.mcmc(post), lag.max = 30)

effectiveSize(as.mcmc(post))

burn.in = 2000
thin_interval = 16
thin_indx = seq(from = burn.in, to = length(post[,1]), by = thin_interval)

thin.post = post[thin_indx, ]
colnames(thin.post) = c("beta1", "beta2", "beta3", "beta4", "beta5", 
                   "beta6", "beta7", "beta8", "beta9", "beta10")

trace.plot(t(thin.post[,1]))
plot(as.mcmc(thin.post))
# coda::autocorr.plot(as.mcmc(thin.post), lag.max = 30)
effectiveSize(as.mcmc(thin.post))
length(thin_indx)

tail(thin.post)



