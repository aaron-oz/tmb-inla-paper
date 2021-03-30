## ######
## ######
## TMB ##
## ######
## ######

message('---- ON SCRIPT 4: running TMB')

## update the tracker
write.table(x=matrix(c(sim.loop.ct, 4), ncol=2), append=T,
            file = paste0(jobtrack.dir,
                          sprintf('exp_%06d_iter_%06d.csv', exp.lvid, exp.iter)), sep=',',
            row.names=F, col.names = F)

## ########
## SETUP ##
## ########

## build design mats for int and covs
## this is done seperately to allow indep. turning each on/off
X_alpha <- matrix(rep(1, nrow(dt), ncol = 1))
X_betas <- as.matrix(dt[, covs[, name], with=FALSE])

templ <- "model_space"
setwd("/homes/azimmer/tmb_inla_comp/")
## TMB::compile(paste0('./', templ,".cpp"))
dyn.load( dynlib(templ) )


## function to convert from data lik string to integer
## allows easily adding more options even though overkill for just 2
tmb.lik.dict <- function(x){
  dict <- list(normal = 0,
               binom = 1)
  dict[[x]]
}

## setup data to feed into the model
data_full <- list(num_i = nrow(dt),  # Total number of observations
                  num_s = mesh_s$n,  # Number of vertices in SPDE mesh
                  y_i   = dt[,Y],    # Number of observed deaths in the cluster
                  n_i   = dt[,N],    # Number of exposures in the cluster
                  X_alpha  = X_alpha,# Covariate design matrix
                  X_betas  = X_betas,# Covariate design matrix
                  M0    = spde$param.inla$M0, # SPDE sparse matrix
                  M1    = spde$param.inla$M1, # SPDE sparse matrix
                  M2    = spde$param.inla$M2, # SPDE sparse matrix
                  Aproj = A.proj,             # Projection matrix
                  options = c(1, ## if 1, run adreport
                              1, ## if 1, use priors
                              ifelse(is.null(alpha), 0, 1),     ## if 1, run with intercept
                              ifelse(is.null(betas), 0, 1),     ## if 1, run with covs
                              ifelse(is.null(clust.var), 0, 1), ## if 1, run with cluster
                              tmb.lik.dict(data.lik),           ## if 0, normal data. if 1, binom data lik
                              1 ## use normalization trick?
                              ),
                  flag = 1, # normalization flag. when 0, prior is returned. when 1 data is included in jnll
                  norm_prec_pri = norm.prec.pri, ## gamma on log(prec)
                  clust_prec_pri = clust.prec.pri, ## gamma on log(prec)
                  alphaj_pri = alphaj.pri, ## normal
                  ## logtau_pri = spde.theta1.pri, ## logtau: normal(mean, prec)
                  ## logkappa_pri = spde.theta1.pri ## logkappa: normal(mean, prec)
                  matern_pri = matern.pri ## c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
                  )

## Specify starting values for TMB parameters for GP
tmb_params <- list(alpha = 0.0, # intercept
                   betas = rep(0, ncol(X_betas)), # cov effects
                   ## log_gauss_sigma = -2, # log(data sd) if using normal dist
                   log_gauss_prec  = -2*-2, # log(data sd)*-2 to get log(prec) if using normal dist
                   log_tau   = -1.0, # Log inverse of tau (Epsilon)
                   log_kappa = -1.0, # Matern range parameter
                   ## log_clust_sigma = -2, # log of cluster sd
                   log_clust_prec = -2*-2, # (log of cluster sd)*-2 to get log(prec)
                   clust_i = rep(0, nrow(dt)), # vector of cluster random effects
                   Epsilon_s = matrix(0, nrow=nodes, ncol=1) # GP value at obs locs
                   )

## make a list of things that are random effects
rand_effs <- c('Epsilon_s')

## NULL out params that aren't in this run and identify extra rand effs if using them
ADmap <- list()

## intercept
if(is.null(alpha)){
  ADmap[['alpha']] <- factor(NA)
}

## fixed effects
if(is.null(betas)){
  ADmap[['betas']] <- rep(factor(NA), ncol(X_betas))
}

## cluster RE
if(is.null(clust.var)){
  ADmap[['log_clust_prec']] <- factor(NA)
  ADmap[['clust_i']] <- rep(factor(NA), nrow(dt))
}else{
  rand_effs <- c(rand_effs, 'clust_i')
}

## normal data obs variance
if(data.lik == 'binom'){
  ADmap[['log_gauss_prec']] <- factor(NA)
}

## make the autodiff generated liklihood func & gradient
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 map = ADmap,
                 hessian=TRUE,
                 DLL=templ)

## should we use the normalization flag?
if(data_full$options[7] == 1){
  obj <- normalize(obj, flag="flag", value = 0) ## value: Value of 'flag' that signifies to not include the data term.
}

## Run optimizer
message('------ fitting TMB')
ptm <- proc.time()[3]

opt0 <- try(do.call("nlminb",
                    list(start       =    obj$par,
                         objective   =    obj$fn,
                         gradient    =    obj$gr,
                         lower = rep(-10, length(obj$par)), ## TODO
                         upper = rep( 10, length(obj$par)), ## TODO
                         control     =    list(trace=1))))
fit_time_tmb <- proc.time()[3] - ptm

if(class(opt0) == "try-error"){
  ## make a flag to see if tmb converged and also check that it converge with a PD cov structure

  ## did the optimizer work?
  ## if not, stop here and record as failure
  tmb.pd.converge <- FALSE
}else{
  ## we can move onto final fit steps and predict
  ## but we still need to test in a minute (line170) if we have PD cov
  tmb.pd.converge <- TRUE

  ## Get standard errors
  SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                      bias.correct = bias.correct,
                      bias.correct.control = list(sd = sd.correct))
  tmb_total_fit_time <- proc.time()[3] - ptm
  tmb_sdreport_time <-  tmb_total_fit_time - fit_time_tmb
  summary(SD0, 'report')

  ## ##########
  ## PREDICT ##
  ## ##########
  message('------ making TMB predictions')

  ## now we can take draws and project to space-time raster locs
  mu <- c(SD0$par.fixed,SD0$par.random)

  ## simulate draws
  ptm2 <- proc.time()[3]
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }

  L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)
  if(class(L) != "try-error"){
    ## then we have a PD precision and we're good to go

    tmb_draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = ndraws)
    tmb_get_draws_time <- proc.time()[3] - ptm2

    ## separate out the tmb_draws
    parnames <- c(names(SD0$par.fixed), names(SD0$par.random))
    epsilon_tmb_draws  <- tmb_draws[parnames == 'Epsilon_s',]
    alpha_tmb_draws    <- matrix(tmb_draws[parnames == 'alpha',], nrow = 1)
    betas_tmb_draws    <- tmb_draws[parnames == 'betas',]
    if(!is.matrix(betas_tmb_draws)) betas_tmb_draws <- matrix(betas_tmb_draws, nrow = 1)
    log_kappa_tmb_draws <- tmb_draws[parnames == 'log_kappa',]
    log_tau_tmb_draws  <- tmb_draws[parnames == 'log_tau',]
    ## log_clust_sigma_tmb_draws <- tmb_draws[parnames == 'log_clust_sigma', ]
    ## log_gauss_sigma_tmb_draws <- tmb_draws[parnames == 'log_gauss_sigma', ]
    log_clust_prec_tmb_draws <- tmb_draws[parnames == 'log_clust_prec', ]
    log_gauss_prec_tmb_draws <- tmb_draws[parnames == 'log_gauss_prec', ]
    sp_range_tmb_draws <- sqrt(8) / exp(log_kappa_tmb_draws)
    sp_sigma_tmb_draws <- (4 * pi * exp(2*log_kappa_tmb_draws) * exp(2 * log_tau_tmb_draws)) ^ (-.5)

    ## values of S at each cell (long by nperiods)
    ## rows: pixels, cols: posterior draws
    pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)

    ## is we have an intercept and no betas, add intercept here
    if(!is.null(alpha) & is.null(betas)){
      ## add on intercept, one alpha draw per row
      pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')
    }

    ## add betas (and intercept if applicable)
    if(!is.null(betas)){

      ## add column for intercept if included
      if(!is.null(alpha)){
        betas_tmb_draws <- rbind(alpha_tmb_draws, betas_tmb_draws)
      }

      ## add on covariate values by draw
      tmb_vals <- list()
      for(p in 1:nperiods) tmb_vals[[p]] <- cov_vals[[p]] %*% betas_tmb_draws

      cell_b <- do.call(rbind, tmb_vals)

      ## add together linear and st components
      pred_tmb <- cell_b + pred_tmb
    }

    ## save prediction timing
    totalpredict_time_tmb <- proc.time()[3] - ptm2

    ## #######
    ## SAVE ##
    ## #######

    # save the fitted model object
    saveRDS(file=sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_fitted_tmb_model.rds',
                         out.dir, exp.lvid, exp.iter),
            object = SD0)

    ## save the posterior param draws
    # saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_param_draws.rds',
    #                        out.dir, exp.lvid, exp.iter),
    #         object = tmb_draws)

    ## save the cell preds
    # saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds.rds',
    #                        out.dir, exp.lvid, exp.iter),
    #         object = pred_tmb)

    ## summarize the latent field
    summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                      sd     = (apply(pred_tmb, 1, sd)))

    ras_med_tmb <- insertRaster(simple_raster, matrix(summ_tmb[, 1], ncol = nperiods))
    ras_sdv_tmb <- insertRaster(simple_raster, matrix(summ_tmb[, 2], ncol = nperiods))

    writeRaster(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds_median_raster.tif',
                               out.dir, exp.lvid, exp.iter),
                x = ras_med_tmb, format='GTiff', overwrite = TRUE)
    writeRaster(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds_stdev_raster.tif',
                               out.dir, exp.lvid, exp.iter),
                x = ras_sdv_tmb, format='GTiff', overwrite = TRUE)

    if(data.lik == 'binom'){
      ## convert to prevalence space and summarize, rasterize, and save again
      pred_tmb_p <- plogis(pred_tmb)

      # saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds_PREV.rds',
      #                        out.dir, exp.lvid, exp.iter),
      #         object = pred_tmb_p)

      summ_tmb_p <- cbind(median = (apply(pred_tmb_p, 1, median)),
                          sd     = (apply(pred_tmb_p, 1, sd)))

      ras_med_tmb_p <- insertRaster(simple_raster, matrix(summ_tmb_p[, 1], ncol = nperiods))
      ras_sdv_tmb_p <- insertRaster(simple_raster, matrix(summ_tmb_p[, 2], ncol = nperiods))

      writeRaster(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds_median_raster_PREV.rds',
                                 out.dir, exp.lvid, exp.iter),
                  x = ras_med_tmb_p, format='GTiff', overwrite = TRUE)
      writeRaster(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_tmb_preds_stdev_raster_PREV.rds',
                                 out.dir, exp.lvid, exp.iter),
                  x = ras_sdv_tmb_p, format='GTiff', overwrite = TRUE)
    }
  }else{
    tmb.pd.converge <- FALSE ## the jointPrec was not PD
    message('------ WARNING: TMB PRECISION IS NOT! PD - mapping to nearest PD precision ')
    message('------ WARNING: THIS RUN WILL BE LOGGED AS TMB FAILING TO CONVERGE ')
  }
}
