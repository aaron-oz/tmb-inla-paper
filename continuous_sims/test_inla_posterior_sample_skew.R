# from geir-arne

set.seed(1020)
x = rnorm(5, mean = 0)
y = rpois(5, exp(x))
res = inla(y~x,
           data = list(x  = x,
                       y = y),
           control.inla = list(int.strategy = "eb"),
           control.compute = list (config = TRUE),
           family = "poisson")
samp = inla.posterior.sample(n = 100000, res, use.improved.mean = T, skew.corr = T)
beta = rep(0, 100000)
for(i in 1:100000)
  beta[i] = samp[[i]]$latent[res$misc$configs$contents$start[3]]
hist(beta, 40, freq = F)
lines(res$marginals.fixed$x)
