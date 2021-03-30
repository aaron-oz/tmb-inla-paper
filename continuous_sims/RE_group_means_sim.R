## this file simulates data from a simple group means RE model and
## calculates the coverage of the group RE for different nominal
## coverages

## the goal is to see if they, on average, capture the coverage

## aoz Sat Jan 18 15:34:33 PST 2020

##############
## libraries #
##############
require(data.table)
require(INLA)
inla.setOption("pardiso.license", "~/sys/licenses/pardiso.lic")
inla.pardiso.check()
require(ggplot2)

sim.data <- function(n.grp = 10,
                     re.s = 1, ## re sd
                     ob.s = 1, ## obs sd if normal liklihood
                     ob.m = 35, ## mean number of observations
                     int = -1, ## intercept
                     lik = 'normal' ## or binom
                     ){

  ## make the random effects of the groups
  re <- rnorm(n = n.grp, mean = 0, sd = re.s)

  ## sample the data
  d <- data.table(grp = integer(),
                  obs = numeric())

  for(i in 1:n.grp){

    n.obs <- rpois(n = 1, lambda = ob.m)

    if(lik == 'normal'){
      obs <- rnorm(n.obs, mean = (int + re[i]), sd = re.s)
    }
    if(lik == 'binomial'){
      obs <- rbinom(n = n.obs, size = 1, p = plogis(int + re[i]))
    }
    
    d <- rbind(d,
               data.table(grp     = i,
                          obs     = obs))

    }
    return(list(data = d,
                true.re = re))
  
}

fit.make.draws <- function(data, ## data from sim.data
                           lik = 'normal',
                           n.draws = 500){

  if(lik == 'normal'){
    mod.fit <- try(suppressWarnings(inla(obs ~ 1 +
                    f(grp, model = 'iid',
                      hyper = list(prec = list(prior = 'pc.prec',
                                               params = c(5, 0.01)))),
                    data = data, family = 'normal',
                    control.compute = list(config = TRUE),
                    control.predictor = list(compute = FALSE))))
  }
  if(lik == 'binomial'){
    mod.fit <- try(suppressWarnings(inla(obs ~ 1 +
                    f(grp, model = 'iid',
                      hyper = list(prec = list(prior = 'pc.prec',
                                               params = c(5, 0.01)))),
                    data = data, family = 'binomial',
                    Ntrials = 1, 
                    control.compute = list(config = TRUE),
                    control.predictor = list(compute = FALSE))))
  }

  ## if INLA errors for any reason, return an 'Error' string
  if(class(mod.fit) == "try-error"){
    draws <- list('Error')
    return('Error')
  }
     
  ## Otherwise, get and return RE draws
  draws <- inla.posterior.sample(n = n.draws, result = mod.fit)
  grp.idx <- grep('grp', rownames(draws[[1]]$latent))
  grp.draws <- do.call('cbind', lapply(draws, function(x){x$latent[grp.idx]}))

  ## return a matrix, long on group RE estimates, wide on draws
  grp.draws 
}

check.coverage <- function(draws,
                           coverage_probs = c(25, 50, 80, 90, 95),
                           true.re){

  d <- data.table(truth = true.re)

  ## get the lower and upper quantiles for each RE for different coverages
  lui <- apply(draws, 1, quantile, 
               p = c((1 - coverage_probs/100)/2, 
                     coverage_probs/100 + (1 - coverage_probs/100) / 2), na.rm=T)
  for(cc in 1:length(coverage_probs)){
    c <- coverage_probs[cc]
    ## message(paste0('-------- calcing ', c ,'% coverage of binomial prob'))
    li       <- lui[cc,]
    ui       <- lui[length(coverage_probs) + cc, ]
    d[, paste0('re.cov.',c)] <- d[['truth']] >= li & d[['truth']] <= ui
  }

  ## average coverage across all REs
  avg.cov <- data.table(true.cov = coverage_probs,
                        obs.cov  = d[,colMeans(.SD), .SDcols = grep('re.cov', names(d))],
                        n.grp    = length(true.re),
                        n.draws  = ncol(draws))
  
}

avg.cov <- data.table(true.cov = integer(),
                      obs.cov  = numeric(), 
                      n.grp    = integer(), 
                      n.draws  = integer(),
                      lik      = character(),
                      num.obs = integer())
for(lik in c('normal', 'binomial')){
  message(sprintf('on lik: %s', lik))
  for(num.re in c(5, 10, 25, 50)){
    message(sprintf('--on num.re: %i', num.re))
    for(num.obs in c(25, 100, 15 ^ 2)){
      message(sprintf('----on num.obs: %i', num.obs))
      for(i in 1:100){
        if(i%%10 == 0)message(sprintf('------on iter: %i', i))

        ## loop over sim'ing and fit'ing until INLA does NOT error
        draws <- list('Error')
        while(draws == 'Error'){
          data    <- sim.data(n.grp = num.re, lik = lik, ob.m = num.obs)
          true.re <- data$true.re
          data    <- data$data
          draws   <- fit.make.draws(data = data, lik = lik)
        }
        ## if draws[[1]] != 'Error'
        avg.cov <- rbind(avg.cov, cbind(check.coverage(draws = draws,
                                                       true.re = true.re),
                                        lik, num.obs))
      }
    }
  }
}

fwrite(avg.cov, '~/tmb_inla_sim/RE_groups_coverage_sim.csv')
avg.cov <- fread('~/tmb_inla_sim/RE_groups_coverage_sim.csv')

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.80, .drop=TRUE) {
  library(plyr)
  
  ## New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     med  = median (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     l.ci = stats::quantile(xx[[col]], probs = (1-conf.interval)/2, na.rm=na.rm), 
                     u.ci = stats::quantile(xx[[col]], probs = 1-(1-conf.interval)/2, na.rm=na.rm) 
                   )
                 },
                 measurevar
  )
  
  ## Rename the "mean" column    
  old.names <- c("mean", 
                 paste0('l.ci.', (1-conf.interval)/2*100, '%'), 
                 paste0('u.ci.', (1-(1-conf.interval)/2)*100, '%'))
  new.names <- c(measurevar, 'l.ci', 'u.ci')
  names(datac)[match(old.names, names(datac))] <- new.names
  
  ## get with width of the ci
  datac$w.ci <- datac$u.ci - datac$l.ci
  
  return(datac)
}

cov.width <- .5
avg.cov.sum <-  summarySE(avg.cov, measurevar="obs.cov",
                           groupvars = c('true.cov', 'n.grp', 'lik', 'num.obs'), 
                           conf.interval = cov.width)
## groups to plot lines
avg.cov.sum$line_group <- paste(avg.cov.sum$lik, avg.cov.sum$n.grp, sep='_')

pd <- position_dodge(0.01)
fit_coverage_CI_summary <- ggplot(avg.cov.sum,
                                  aes(x = true.cov / 100, y = obs.cov,
                                      shape = lik, 
                                      linetype = lik,
                                      color = as.factor(n.grp),
                                      group = line_group),
                                  position=position_jitter(w=0.02, h=0.02)) + 
geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
  facet_wrap(~ num.obs, scales = "fixed") +
geom_line(position=pd) +
geom_point(position=pd) +
geom_abline(intercept = 0, slope = 1) +
ggtitle('sim coverage from RE groups') +
## fix the legends a bit
labs(color = "num.grps", shape='lik', linetype='lik') + ## legend titles
xlab('Nominal Coverage') + 
ylab('Observed Monte Carlo Coverage')

fit_coverage_CI_summary

ggsave('~/GeneralExam/RE_grp_cov_sim.png', 
        plot = fit_coverage_CI_summary ,
        device = 'png', units = 'in',
        width = 12, height = 8)
