## #######
## #######
## INLA ##
## #######
## #######

message('---- ON SCRIPT 5: running INLA')

## update the tracker
write.table(x=matrix(c(sim.loop.ct, 5), ncol=2), append=T,
            file = paste0(jobtrack.dir,
                          sprintf('exp_%06d_iter_%06d.csv', exp.lvid, exp.iter)), sep=',',
            row.names=F, col.names = F)

## ########
## SETUP ##
## ########

A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = dt.coords,
  group = dt.pers
)
space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde,
                                n.group = nperiods)

inla.covs <- covs$name
design_matrix <- data.frame(int = rep(1, nrow(dt)), dt[, inla.covs, with=F], clust.id = 1:nrow(dt))
stack.obs <- inla.stack(tag='est',
                        data=list(Y=dt$Y), ## response
                        A=list(A,1), ## proj matrix, not sure what the 1 is for
                        effects=list(
                          space,
                          design_matrix)
                        )

formula <- formula(paste('Y ~ -1',
                         ifelse(is.null(alpha), '', ' + int'),
                         ifelse(is.null(betas), '', paste0(' + ', (paste(inla.covs, collapse = ' + ')))),
                         ifelse(is.null(clust.var), '',
                                paste0(' + f(clust.id, model = \'iid\', hyper = list(prec = list(prior=\'pc.prec\', param = c(',
                                       clust.prec.pri[1],', ', clust.prec.pri[2], '))))')),
                         ' + f(space, model = spde, group = space.group, control.group = list(model = \'ar1\'))',
                         sep = ''))

## function to convert from data lik string to integer
## allows easily adding more options even though overkill for just 2
inla.lik.dict <- function(x){
  dict <- list(normal = 'normal',
               binom = 'binomial')
  dict[[x]]
}
inla.setOption("enable.inla.argument.weights", TRUE)



## ######
## FIT ##
## ######
message('------ fitting INLA')

ptm <- proc.time()[3]
if(data.lik == 'normal'){
  res_fit <- try(inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       mean = list(default = alphaj.pri[1]),          ## fixed effects prior mean
                                       prec = list(default = 1 / alphaj.pri[2] ^ 2)), ## fixed effects prior prec
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat),
                  control.compute=list(config = TRUE),
                  control.family = list(hyper = list(prec = list(prior="pc.prec", param = norm.prec.pri))),
                  family = inla.lik.dict(data.lik),
                  num.threads = cores, #
                  scale = dt$N,
                  verbose = FALSE, ## this must be false to get the logfile
                  keep = FALSE))
}else if(data.lik == 'binom'){
  res_fit <- try(inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           ## link = 1, ## removed after looking at NMM
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec = list(default = 1 / alphaj.pri[2] ^ 2)),
                  control.inla = list(strategy = inla.approx,
                                      int.strategy = inla.int.strat ##,
                                      ## h = 1e-3, ## removed after looking at NMM
                                      ## tolerance = 1e-6 ## removed after looking at NMM
                                      ),
                  control.compute=list(config = TRUE),
                  family = inla.lik.dict(data.lik),
                  num.threads = cores, #
                  Ntrials = dt$N,
                  verbose = FALSE, ## this must be false to get the logfile!
                  keep = FALSE))
}
fit_time_inla <- proc.time()[3] - ptm

## check to see if INLA crashed
if(class(res_fit) == "try-error"){
  inla.crash <- TRUE
  inla.pd.hess <- FALSE       ## set this since otherwise since it's assignment will be skipped if INLA crashes
  inla.mode.converge <- FALSE ## set this since otherwise since it's assignment will be skipped if INLA crashes
}else{
  inla.crash <- FALSE
}

if(!inla.crash){

  ## and grep the logfile for the hessian pd error where INLA self-corrects by adding to the diagonal
  ## the full warning looks something like this:
  # *** WARNING *** Eigenvalue 2 of the Hessian is -0.183833 < 0
  # *** WARNING *** Set this eigenvalue to 2.03328
  # *** WARNING *** This have consequence for the accurancy of
  # *** WARNING *** the approximations; please check!!!
  # *** WARNING *** R-inla: Use option inla(..., control.inla = list(h = h.value), ...)
  # *** WARNING *** R-inla: to chose a different `h.value'.
  inla.pd.hess <- ifelse(sum(grepl('Eigenvalue . of the Hessian', res_fit$logfile)) == 0, TRUE, FALSE)

  ## check to see if INLA converged nicely
  ## using the check suggested here: https://groups.google.com/forum/#!topic/r-inla-discussion-group/LsCpuCsr-Qo
  ## and noting that failing this check may not be terrible
  inla.mode.converge <- ifelse(res_fit$mode$mode.status == 0, TRUE, FALSE)

  ## if the model didn't crash, AND it the hess was PD, AND if the mode converged,
  ## then we predict out
  ## we could predict if either of these failed, but I'm considering these convergence failures
  if(TRUE){# inla.pd.hess & inla.mode.converge){ ## added Mar 14, 2020

    ## ##########
    ## PREDICT ##
    ## ##########
    message('------ making INLA predictions')

    ptm <- proc.time()[3]
    inla_draws <- inla.posterior.sample(ndraws, res_fit, use.improved.mean = bias.correct)
    inla_get_draws_time <- proc.time()[3] - ptm

    ## get parameter names
    par_names <- rownames(inla_draws[[1]]$latent)

    ## index to spatial field and linear coefficient samples
    s_idx <- grep('^space.*', par_names)
    l_idx <- which(!c(1:length(par_names)) %in% grep('^space.*|Predictor|clust.id', par_names))
    ## deleted |[*:*] on 14Mar2020 from grep call

    ## get spatial draws as matrices and project to deaws at locations
    pred_s <- sapply(inla_draws, function (x) x$latent[s_idx])
    pred_inla <- as.matrix(A.pred %*% pred_s)

    ## get intercept and coef draws and convert to covariate effects
    if(!is.null(alpha) | !is.null(betas)){
      pred_l <- sapply(inla_draws, function (x) x$latent[l_idx])
      if(!is.matrix(pred_l)){
        pred_l <- matrix(pred_l, ncol = length(pred_l))
      }
      rownames(pred_l) <- res_fit$names.fixed

      ## extract cell values  from covariates, deal with timevarying covariates here
      non_space_names <- par_names[l_idx]
      ## remove :N if present
      non_space_names <- unlist(lapply(non_space_names,
                                function(x){
                                  if(grepl(':', x)){
                                    x <- strsplit(x, ':')[[1]][1]
                                  }
                                  return(x)
                                }))
      cov_effs <- list()
      for(p in 1:nperiods)  cov_effs[[p]] <- cov_vals[[p]][, non_space_names] %*% pred_l

      cov_effs <- do.call(rbind, cov_effs)

      pred_inla <- pred_inla + cov_effs
    }

    if(!is.null(clust.var)){
      ## get draws of clust precision
      clust_prec_inla_draws <- sapply(inla_draws, function(x) {
        clust.idx <- which(grepl('clust.id', names(inla_draws[[1]]$hyper)))
        x$hyperpar[[clust.idx]]}) ## this gets the precision for the cluster RE
    }

    totalpredict_time_inla <- proc.time()[3] - ptm

    ## #######
    ## SAVE ##
    ## #######

    ## save the fitted model object
    saveRDS(file=sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_fitted_inla_model.rds',
                                                out.dir, exp.lvid, exp.iter),
            object = res_fit)

    ## save the posterior param draws
    ## saveRDS(file = sprintf('%s/modeling/outputs/tmb/experiment%06d_iter%06d_inla_param_draws.rds',
    ##                       out.dir, exp.lvid, exp.iter),
    ##         object = inla_draws)

    ## save the cell_pred
    # saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds.rds',
    #                        out.dir, exp.lvid, exp.iter),
    #         object = pred_inla)

    ## summarize the latent field
    summ_inla <- cbind(median = (apply(pred_inla, 1, median)),
                       sd     = (apply(pred_inla, 1, sd)))

    ras_med_inla <- insertRaster(simple_raster, matrix(summ_inla[, 1], ncol = nperiods))
    ras_sdv_inla <- insertRaster(simple_raster, matrix(summ_inla[, 2], ncol = nperiods))

    writeRaster(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds_median_raster.tif',
                               out.dir, exp.lvid, exp.iter),
                x = ras_med_inla, format='GTiff', overwrite=TRUE)
    writeRaster(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds_stdev_raster.tif',
                               out.dir, exp.lvid, exp.iter),
                x = ras_sdv_inla, format='GTiff', overwrite=TRUE)

    if(data.lik == 'binom'){
      ## convert to prevalence space and summarize, rasterize, and save again
      pred_inla_p <- plogis(pred_inla)

      # saveRDS(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds_PREV.rds',
      #                        out.dir, exp.lvid, exp.iter),
      #         object = pred_inla_p)

      summ_inla_p <- cbind(median = (apply(pred_inla_p, 1, median)),
                           sd     = (apply(pred_inla_p, 1, sd)))

      ras_med_inla_p <- insertRaster(simple_raster, matrix(summ_inla_p[, 1], ncol = nperiods))
      ras_sdv_inla_p <- insertRaster(simple_raster, matrix(summ_inla_p[, 2], ncol = nperiods))

      writeRaster(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds_median_raster_PREV.rds',
                                 out.dir, exp.lvid, exp.iter),
                  x = ras_med_inla_p, format='GTiff', overwrite=TRUE)
      writeRaster(file = sprintf('%s/modeling/outputs/inla/experiment%06d_iter%06d_inla_preds_stdev_raster_PREV.rds',
                                 out.dir, exp.lvid, exp.iter),
                  x = ras_sdv_inla_p, format='GTiff', overwrite=TRUE)
    }
  }## inla.pd.hess & inla.mode.converge
} ## !inla.crash
