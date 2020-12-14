rm(list=ls())

library(coda)
library(R.matlab)
library(spate)
setwd("/Users/stevehof/Documents/School/Fall_2020/STAT_460/Assignments/Quality_Work")


dat = readMat("Dataset_FinalProject_2.mat")

X = dat$X
y = dat$Y

indx_gen = function(j, M) {
  n = dim(M)[2]
  if (j == 1) return(M[1, ])
  if (j == n) return(M[2, ])
  return(c(M[2, 1 : j-1], M[1, j:n]))
}

update_pi = function(delta_j, a_pi, b_pi) {
  shape1 = sum(delta_j) + a_pi/2
  shape2 = sum(1-delta_j) + b_pi/2
  rbeta(n = 1, shape1, shape2)
}

update_delta = function(bet, pi, tau, eps, j) {
  p0 = ((1 - pi) * tau / sqrt(eps)) * exp(-1/(2*eps) * bet[j]^2)
  p1 = pi * exp(-1/(2*tau^2) * bet[j]^2)
  rob = p1 / (p0 + p1)
  return(rbinom(n = 1, size = 1, prob = rob))
  
}

update_beta = function(delta_j, j, bet, tau, eps, mutop = 0, mubot = 0) {
  partial.mutop = mutop
  partial.mubot = mubot
  for(i in 1:dim(X)[1]) {
    partial.mutop = partial.mutop + X[i, j] * (y[i] - sum(setdiff(X[i,], X[i,j]) * 
                                                            setdiff(bet, bet[j])))
    partial.mubot = partial.mubot + X[i, j]^2
  }
  
  if(delta_j == 0) {
    bot = 1 + eps * partial.mubot
    mu = eps * partial.mutop / bot
    sd = sqrt(eps * bot^-1)
    return(rnorm(n = 1, mean = mu, sd = sd))
    
  }
  if(delta_j == 1) {
    bot = 1 + tau^2 * partial.mubot
    mu = tau^2 * partial.mutop / bot
    sd = sqrt(tau^2 * bot^-1)
    return(rnorm(n = 1, mean = mu, sd = sd))
  }
}

gibbs = function(n_iter, init, priors) {
  # initialize
  beta.out = matrix(data = NA, nrow = n_iter, ncol = dim(X)[2])
  delta.curr = init$delta_j
  beta.curr = init$beta
  beta.out[1, ] = beta.curr
  
  # Gibbs sampler
  for(k in 2:n_iter) {
    pi.curr = update_pi(delta_j = delta.curr, a_pi = priors$a_pi, 
                        b_pi = priors$b_pi)
    
    for(j in 1:length(beta.curr)) {
      
      delta.curr = update_delta(bet = beta.out[k-1,], pi = pi.curr, 
                                tau = priors$tau, eps = priors$eps, j)

      betas = indx_gen(j = j, M = beta.out[(k-1) : k, ])
      beta.curr[j] = update_beta(delta_j = delta.curr, j = j, bet = betas, 
                                 tau = priors$tau, eps = priors$eps)
      
      beta.out[k, j] = beta.curr[j]
    }
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

priors$a_pi = 12
priors$b_pi = 1
priors$tau = 10
priors$eps = 10^-4
# postone12one = gibbs(n_iter = n_iter, init = init, priors = priors)
# post = readRDS(file = "Quality_Work/post_1.rds")

# colnames(postone12one) = c("beta1", "beta2", "beta3", "beta4", "beta5", 
#                    "beta6", "beta7", "beta8", "beta9", "beta10")
# saveRDS(postone12one, file = "postone12one.Rds")

post1450 = readRDS(file = "postone1450.Rds")
post12 = readRDS(file = "post_1_2.Rds")
postone12 = readRDS(file = "postone12.Rds")
post12one = readRDS(file = "postone12one.Rds")
all.posts.new.pi = list(post1450, post12, postone12, post12one)

burn.and.thin = function(post, bi, ti) {
  thin.indx = seq(from = bi, to = length(post[, 1]), by = ti)
  thin.post = post[thin.indx, ]
  colnames(thin.post) = c("beta1", "beta2", "beta3", "beta4", "beta5",
                          "beta6", "beta7", "beta8", "beta9", "beta10")
  return(thin.post)
}
thin.all.post.new = list(4)
for (i in 1:4) {
  thin.all.post.new[[i]] = burn.and.thin(post = all.posts.new.pi[[i]], bi = 2000, ti = 10)
}
thin.post1450 = thin.all.post.new[[1]]
thin.post12 = thin.all.post.new[[2]]
thin.postone12 = thin.all.post.new[[3]]
thin.post12one = thin.all.post.new[[4]]
# tail(post)
# plot(as.mcmc(post))
bob = summary(as.mcmc(thin.post1450))
bob$quantiles
bob$statistics

barb = t(hdi(as.mcmc(thin.post1450)))
barb
# coda::autocorr.plot(as.mcmc(post), lag.max = 30)
# 
# effectiveSize(as.mcmc(post))
# 
# burn.in = 2000
# thin_interval = 20
# thin_indx = seq(from = burn.in, to = length(post[,1]), by = thin_interval)
# 
# thin.post = post[thin_indx, ]
# colnames(thin.post) = c("beta1", "beta2", "beta3", "beta4", "beta5", 
#                         "beta6", "beta7", "beta8", "beta9", "beta10")


plot(as.mcmc(thin.post))
# coda::autocorr.plot(as.mcmc(thin.post), lag.max = 30)
effectiveSize(as.mcmc(thin.post))
length(thin_indx)
coda::autocorr.plot(as.mcmc(thin.post), lag.max = 30)




