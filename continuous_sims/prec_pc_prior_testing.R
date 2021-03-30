## need to try out a bunch of pc.priors...


prior.us <- c(.1, .25, .5, .75, 1)
prior.as <- c(.01, .05, .1, .25, .5)

pc.test.res <- data.table(value = numeric(),
                          param = character(),
                          model = character(),
                          tmb.c = logical(),
                          u     = numeric(),
                          a     = numeric())

for(uu in 1:length(prior.us)){
  for(aa in 1:length(prior.as)){
    prior.u <- prior.us[uu]
    prior.a <- prior.as[aa]
    
    norm.prec.pri <- clust.prec.pri <- c(prior.u, prior.a)
    
    source('./4_setup_run_predict_tmb.R')
    
    ## update convergence args for while loop
    if(tmb.pd.converge){
      tmb.converge <- 1
    }else{
      tmb.converge <- 0
      tmb.converge.fails <- tmb.converge.fails + 1
    }
    
    ## #######
    ## INLA ##
    ## #######
    source('./5_setup_run_predict_inla.R')
    
    ## update convergence args for while loop
    if(inla.mode.converge){
      inla.converge <- 1
    }else{
      inla.converge <- 0
      inla.converge.fails <- inla.converge.fails + 1
    }
    
    pc.test.res <- rbind(pc.test.res, data.table(value = unname(SD0$value['clust_prec']),
                                                 param = 'clust_pres',
                                                 model = 'TMB',
                                                 tmb.c = as.logical(tmb.converge),
                                                 u     = prior.u,
                                                 a     = prior.a))
    
    pc.test.res <- rbind(pc.test.res, data.table(value = unname(SD0$value['gauss_prec']),
                                                 param = 'gauss_pres',
                                                 model = 'TMB',
                                                 tmb.c = as.logical(tmb.converge),
                                                 u     = prior.u,
                                                 a     = prior.a))
    
    pc.test.res <- rbind(pc.test.res, data.table(value = res_fit$summary.hyperpar[2, 4],
                                                 param = 'clust_pres',
                                                 model = 'INLA',
                                                 tmb.c = as.logical(tmb.converge),
                                                 u     = prior.u,
                                                 a     = prior.a))
    
    pc.test.res <- rbind(pc.test.res, data.table(value = res_fit$summary.hyperpar[1, 4],
                                                 param = 'gauss_pres',
                                                 model = 'INLA',
                                                 tmb.c = as.logical(tmb.converge),
                                                 u     = prior.u,
                                                 a     = prior.a))
  }
}

# ## i forgot to include a and u the first time...
# pc.test.res[, u := rep(prior.us, each=length(prior.as)*4)]
# pc.test.res[, a := rep(rep(prior.as, each=4), length(prior.us))]


## TODO compare this (pc.test.res.prec ) last run to the one that ran overnight (pc.test.res)
## pc.test.res.prec <- pc.test.res ## without log(prec) change in tmb

gg <- ggplot(pc.test.res, aes(x=param, y=value, col=model)) + 
  geom_point(aes(pch=tmb.c, size=5)) + 
  facet_grid(u~a) +
  ylim(c(0,1000)) + 
  geom_abline(intercept = 100, slope = 0)
gg




