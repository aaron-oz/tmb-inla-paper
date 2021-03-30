## #############
## #############
## VALIDATION ##
## #############
## #############

message('---- ON SCRIPT 6: running validation')

## update the tracker
write.table(x=matrix(c(sim.loop.ct, 6), ncol=2), append=T,
            file = paste0(jobtrack.dir, 
                          sprintf('exp_%06d_iter_%06d.csv', exp.lvid, exp.iter)), sep=',',
            row.names=F, col.names = F)

## 1) summarize fitted params
## 2) big plots showing difference in fits
## 3) calcualte and summarize predictive metrics  

## ###################################
## 1) summarize fitted param values ##
## ###################################
message('------ making table of true and estimated params')

## make a dt for comparing results and store some relevant computing and mesh params
res <- data.table(st_mesh_nodes = rep(nrow(epsilon_tmb_draws),2))
res[, cores           := rep(cores,2)]
res[, s_mesh_max_edge := rep(eval(parse(text = mesh_s_params))[2],2)]
res[, s_mesh_cutoff   := rep(eval(parse(text = mesh_s_params))[1],2)]
res[, draws           := c(ndraws,ndraws)]

## time variables
res[,fit_time  := c(fit_time_inla,fit_time_tmb)]
res[,pred_time := c(totalpredict_time_inla,totalpredict_time_tmb)]
res[,pt_tmb_sdreport_time := c(NA,tmb_sdreport_time)]
res[,pt_get_draws_time := c(inla_get_draws_time,tmb_get_draws_time)]

## convergence
res[, convergence := c(inla.converge, tmb.converge)]
## converge attempts
res[, convergence.fails := c(inla.converge.fails, tmb.converge.fails)]

## fe coefficients
if(!is.null(alpha) | !is.null(betas)){
  res[, paste0('fe_',res_fit$names.fixed,'_mean') := 
        as.data.frame(rbind(res_fit$summary.fixed$mean, SD0$par.fixed[1:length(res_fit$names.fixed)]))]
  res[, paste0('fe_',res_fit$names.fixed,'_sd')   := 
        as.data.frame(rbind(res_fit$summary.fixed$sd, sqrt(diag(SD0$cov.fixed))[1:length(res_fit$names.fixed)]))]
}

## cluster var
if(!is.null(clust.var)){
  inla.var <- 1/res_fit$summary.hyperpar[grep('clust.id',rownames(res_fit$summary.hyperpar)),'0.5quant']
  tmb.var  <- 1/SD0$value[grep('clust_prec', names(SD0$value))]
  res[,clust_var_mean := unname(c(inla.var, tmb.var))]
  ## apply delta method. var(f(x)) = var(x)*f'(x)^2 -> var(1/prec)  = var(prec)* (-1/prec^2)^2 = var(prec)/prec^4 = sd(prec)/prec^2
  res[,clust_var_sd := c(res_fit$summary.hyperpar[grep('clust.id',
                                                       rownames(res_fit$summary.hyperpar)), 'sd'] * inla.var^2,
                         SD0$sd[grep('clust_prec', names(SD0$value))] * tmb.var^2) ]
}

## normal data obs var
if(data.lik == 'normal'){
  inla.var <- 1/res_fit$summary.hyperpar[grep('Gaussian',rownames(res_fit$summary.hyperpar)),'0.5quant']
  tmb.var  <- 1/SD0$value[grep('gauss_prec', names(SD0$value))]
  res[,gauss_var_mean := unname(c(inla.var, tmb.var))]
  ## apply delta method. var(f(x)) = var(x)*f'(x)^2 -> var(1/prec)  = var(prec)* (-1/prec^2)^2 = var(prec)/prec^4
  res[,gauss_var_sd := c(res_fit$summary.hyperpar[grep('Gaussian',
                                                       rownames(res_fit$summary.hyperpar)), 'sd'] * inla.var^2,
                             SD0$sd[grep('gauss_prec', names(SD0$value))] * tmb.var^2) ]
}

## hyperparameters
# inla.log.tau <- res_fit$summary.hyperpar[grep('Theta1',rownames(res_fit$summary.hyperpar)),4]
# tmb.log.tau  <- SD0$par.fixed['log_tau']
# res[,matern_logtau_mean := unname(c(inla.log.tau, tmb.log.tau))]
# res[,matern_logtau_sd := c(res_fit$summary.hyperpar[grep('Theta1',rownames(res_fit$summary.hyperpar)),2],
#                              sqrt(SD0$cov.fixed['log_tau','log_tau'])) ]
# 
# inla.log.kappa <- res_fit$summary.hyperpar[grep('Theta2',rownames(res_fit$summary.hyperpar)),4]
# tmb.log.kappa  <- SD0$par.fixed['log_kappa']
# res[,matern_logkappa_mean := unname(c(inla.log.kappa, inla.log.tau)) ]
# res[,matern_logkappa_sd := c(res_fit$summary.hyperpar[grep('Theta2',rownames(res_fit$summary.hyperpar)),2],
#                                sqrt(SD0$cov.fixed['log_kappa','log_kappa'])) ]

inla.sp.range <- res_fit$summary.hyperpar[grep('Range for space',
                                               rownames(res_fit$summary.hyperpar)), '0.5quant']
tmb.sp.range  <- SD0$value['sp_range']
res[,matern_range_mean := unname(c(inla.sp.range, tmb.sp.range))]
res[,matern_range_sd := c(res_fit$summary.hyperpar[grep('Range for space',
                                                        rownames(res_fit$summary.hyperpar)),'sd'],
                          SD0$sd[which(names(SD0$value)=='sp_range')]) ]

inla.sp.sigma <- res_fit$summary.hyperpar[grep('Stdev for space',rownames(res_fit$summary.hyperpar)),'0.5quant']
tmb.sp.sigma  <- SD0$value['sp_sigma']
res[,matern_sigma_mean := unname(c(inla.sp.sigma, tmb.sp.sigma))]
res[,matern_sigma_sd := c(res_fit$summary.hyperpar[grep('Stdev for space',
                                                        rownames(res_fit$summary.hyperpar)),'sd'],
                          SD0$sd[which(names(SD0$value)=='sp_sigma')]) ]

## add extra row to filled with the truth
res <- rbind(lapply(1:ncol(res), function(x){NA}), res)

## slot in the truth and also make a list of all params in the model
params <- NULL
if(!is.null(alpha)){ res[1, fe_int_mean := alpha]; params <- c(params, 'alpha')}
if(!is.null(betas) & is.null(alpha)) { 
  res[1, grep('fe.*mean', colnames(res), value=T) := as.list(betas)]
  params <- c(params, rep('beta', length(betas)))}
if(!is.null(betas) & !is.null(alpha)) {
  res[1, grep('fe.*mean', colnames(res), value=T)[-1] := as.list(betas)]
  params <- c(params, rep('beta', length(betas)))}
if(!is.null(clust.var)) {res[1, clust_prec_mean:= 1 / clust.var]; params <- c(params, 'clust.prec')}
if(data.lik == 'normal') {res[1, gauss_prec_mean:= 1 / norm.var]; params <- c(params, 'gauss.prec')}
res[1, matern_range_mean := sp.range];     params <- c(params, 'sp.range')
res[1, matern_sigma_mean := sqrt(sp.var)]; params <- c(params, 'sp.sigma')


rr <- data.table(item=colnames(res))
rr <- cbind(rr, t(res))
names(rr) <- c('quantity','TRUE', 'R-INLA','TMB')
rr$diff <- rr[,3]-rr[,4]

write.table(x = rr, row.names = FALSE, sep=',', 
          file = sprintf('%s/validation/experiment%06d_iter%06d_param_summary_table.csv', out.dir, exp.lvid, exp.iter))

## we can now plot this table with: grid.table(rr)

## ###########################
## 2) make a bunch of plots ##
## ###########################
## pdf(sprintf('%s/validation/inla_tmb_summary_comparison_plot_%i_new.pdf',out.dir, exp.iter), height=15,width=30)
## TODO? one overall plot? or somehow stitch together later...

## ~~~~~~~~~~~~~~~~~~~~~~~~~
## plot the table of results
## ~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting the summary table')

png(sprintf('%s/validation/experiment%06d_iter%06d_plot_01_summary_table.png', out.dir, exp.lvid, exp.iter),
    height=7, width=9, units = 'in', res = 250)
cols <- names(rr)[2:5]
rr[,(cols) := round(.SD, 3), .SDcols=cols]
grid.table(rr)
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot priors and posteriors
## also get coverage for params
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting priors and posteriors')
message('------- (while also getting coverage for each param)')

## assume:
## 1) alphas are normal
## 2) logtau and logkappa are normal
## 3) log(SD) is pc.prior

## NOTE: also assume that intercept is the first 'beta' and that betas are listed first!

png(sprintf('%s/validation/experiment%06d_iter%06d_plot_02_parameter_densities.png', out.dir, exp.lvid, exp.iter),
    height=9, width=9, units = 'in', res = 250)

## set layout: start with a square, and drop rows is not needed
num.dists <- length(params)
par.dim <- rep(ceiling(sqrt(num.dists)), 2)
while((par.dim[1]-1)*par.dim[2] >= num.dists){
  par.dim[1] <- par.dim[1]-1
} 
par(mfrow = par.dim)

## make the plots
for(ii in 1:num.dists){

  param <- params[ii]
  message(sprintf('-------- plotting prior and post for: %s', param))
  
  ## get prior curves and posterior draws
  if(param == 'alpha'){
    param.name <- names(res_fit$marginals.fixed)[ii]

    if(param.name == 'int'){
      true.val <- alpha
    }else{
      true.val <- betas[ii - 1]
    }
    
    prior.mean <- alphaj.pri[1]
    prior.sd   <- alphaj.pri[2]

    tmb.post.draws  <- alpha_tmb_draws[ii, ]
    inla.post.draws <- pred_l[ii, ]

    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
  }
  
  if(param == 'beta'){

    param.name <- names(res_fit$marginals.fixed)[ii]
    if(param.name == 'int'){
      true.val <- alpha
    }else{
      true.val <- betas[ii - 1]
    }
    
    prior.mean <- alphaj.pri[1]
    prior.sd   <- sqrt(alphaj.pri[2])
    
    tmb.post.draws  <- betas_tmb_draws[ii, ]
    inla.post.draws <- pred_l[ii, ]
    
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
  }
  
  if(param == 'logkappa'){
    
    true.val <- logkappa

    mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
    prior.mean <- mesh.pri$theta.prior.mean[2]
    prior.prec <- mesh.pri$theta.prior.prec[2, 2]
    prior.sd   <- 1 / sqrt(prior.prec)

    tmb.post.draws  <- log_kappa_tmb_draws
    inla.post.dist  <- res_fit$marginals.hyperpar[['Theta2 for space']]
    ## to sample from the approximate posterior, we need to get the discrete probabilities for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    }
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

  
    tmb.post.median  <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log kappa"
  }
  
  if(param == 'logtau'){

    true.val <- logtau

    mesh.pri   <- param2.matern.orig(mesh_s) ## theta1 = logtau, theta2 = logkappa 
    prior.mean <- mesh.pri$theta.prior.mean[1]
    prior.prec <- mesh.pri$theta.prior.prec[1, 1]
    prior.sd   <- 1 / sqrt(prior.prec)
    
    tmb.post.draws <- log_tau_tmb_draws
    inla.post.dist  <- res_fit$marginals.hyperpar[['Theta1 for space']]
    ## to sample from the approximate posterior, we need to get the discrete probabilities for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    }
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- dnorm(x.prior, mean = prior.mean, sd = prior.sd)

    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "log tau"
  }
  
  if(param == 'sp.range'){
    
    true.val <- sp.range
    
    all.hyper <- inla_all_hyper_postprocess(res_fit$all.hyper)
    hyper   <- res_fit$marginals.hyperpar
    hyp.id  <- grep('Range for space', names(hyper))
    id <- unlist(strsplit(attr(hyper[[hyp.id]], "hyperid"), "\\|"))
    prior <- INLA:::inla.extract.prior(tolower(id[2]), id[1], all.hyper)
    l1 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[1]
    l2 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[2]
    d <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[3]
    
    # calculate pc priors as done in the Fuglstad 2017 paper
    calc_pc_prior_ranges <- function(ranges){
      d/2 * l1 * ranges^(-d/2 - 1) * exp(-l1 * ranges^(-d/2))
    }
    
    tmb.post.draws  <- sp_range_tmb_draws
    inla.post.dist  <- res_fit$marginals.hyperpar[['Range for space']]
    ## to sample from the approximate posterior, we need to get the discrete probabilities for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    } 
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- calc_pc_prior_ranges(x.prior)
    
    tmb.post.median  <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "Matern range"
  }
  
  if(param == 'sp.sigma'){
    
    true.val <- sqrt(sp.var)
    
    ## we get the params from the joint pcmatern using the range, even though this is for the sigma/stdev
    all.hyper <- inla_all_hyper_postprocess(res_fit$all.hyper)
    hyper   <- res_fit$marginals.hyperpar
    hyp.id  <- grep('Range for space', names(hyper))
    id <- unlist(strsplit(attr(hyper[[hyp.id]], "hyperid"), "\\|"))
    prior <- INLA:::inla.extract.prior(tolower(id[2]), id[1], all.hyper)
    l1 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[1]
    l2 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[2]
    d <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[3]
    
    calc_pc_prior_sigmas <- function(sigmas){
      l2 * exp(- l2 * sigmas)
    }
    
    tmb.post.draws  <- sp_sigma_tmb_draws
    inla.post.dist  <- res_fit$marginals.hyperpar[['Stdev for space']]
    ## to sample from the approximate posterior, we need to get the discrete probabilities 
    ## for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    }
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior    <- seq(xlim[1], xlim[2], len = 1000)
    y.prior    <- calc_pc_prior_sigmas(x.prior)
    
    tmb.post.median  <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "Matern sigma"
  }
  if(param == 'gauss.prec'){

    true.val <- 1 / norm.var

    prior.u  <- norm.prec.pri[1]
    prior.a <- norm.prec.pri[2]
    
    tmb.post.draws <- exp(log_gauss_prec_tmb_draws)
    inla.post.dist  <- res_fit$marginals.hyperpar[['Precision for the Gaussian observations']]
    ## to sample from the approximate posterior, we need to get the discrete probabilities 
    ## for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    }

    if(true.val==Inf) {
      ## get a safe range for plotting, and get the prior
      xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws),
                              probs=c(.0001, .95)) ## avoid crazy extremes
    }else{
      ## get a safe range for plotting, and get the prior
      xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                              probs=c(.0001, .95)) ## avoid crazy extremes
      xlim <- range(c(xlim, true.val))
    }
    if(xlim[1] == -Inf) xlim[1] <- -1000
    if(xlim[2] == +Inf) xlim[2] <- +1000 
    
    ## this posterior can have crazy long right tail... 
    ## for visualization purposes, we set the max xlim to be 5x the true.val
    if(xlim[2] > 5*true.val){xlim[2] <- 5*true.val}
    
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dPCPriPrec(x.prior, u = prior.u, a=prior.a, give_log = 0)
    
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "Gauss prec"
  }
  if(param == 'clust.prec'){

    true.val <- 1 / clust.var

    prior.u <- clust.prec.pri[1]
    prior.a <- clust.prec.pri[2]

    tmb.post.draws <- exp(log_clust_prec_tmb_draws)
    inla.post.dist  <- res_fit$marginals.hyperpar[[grep('clust.id', names(res_fit$marginals.hyperpar))]]
    ## to sample from the approximate posterior, we need to get the discrete probabilities 
    ## for each x-value. we do this by assuming an approximate rectangle
    inla.post.dist.widths <- diff(inla.post.dist[,1])
    inla.post.dist.widths <- c(inla.post.dist.widths[1] / 2,
                               (inla.post.dist.widths[-length(inla.post.dist.widths)] + inla.post.dist.widths[-1])/2,
                               inla.post.dist.widths[length(inla.post.dist.widths)]/2)
    inla.post.dist.probs <- inla.post.dist.widths * inla.post.dist[,2] ## width * pdf heigh approx= prob
    # drop NaNs, and Infs
    inla.post <- cbind(inla.post.dist[,1], inla.post.dist.widths, inla.post.dist.probs)
    inla.post <- na.omit(inla.post)
    inf.ind   <- unique(c(which(inla.post[,1]==Inf), which(inla.post[,2]==Inf),  which(inla.post[,3]==Inf)))
    if(length(inf.ind) > 0){
      inla.post.dist.x      <- inla.post[-inf.ind, 1]
      inla.post.dist.widths <- inla.post[-inf.ind, 2]
      inla.post.dist.probs  <- inla.post[-inf.ind, 3]
    }else{
      inla.post.dist.x      <- inla.post[, 1]
      inla.post.dist.widths <- inla.post[, 2]
      inla.post.dist.probs  <- inla.post[, 3]
    }
    ## rarely, the inla posterior dist is a point mass (eg pc.prec posterior has point mass at Tau=Inf!)
    if(sum(inla.post.dist.probs)==0){
      inla.post.draws <- sample(x=inla.post.dist[,1], size=ndraws, replace=TRUE, prob=inla.post.dist[,2])
    }else{
      inla.post.draws <- sample(x=inla.post.dist.x, size=ndraws, replace=TRUE, prob=inla.post.dist.probs)
    }
    ## get a safe range for plotting, and get the prior
    xlim <- stats::quantile(c(tmb.post.draws, inla.post.draws, true.val),
                            probs=c(.0001, .9999)) ## avoid crazy extremes
    if(xlim[1] == -Inf) xlim[1] <- -1000; if(xlim[2] == +Inf) xlim[2] <- +1000 
    x.prior <- seq(xlim[1], xlim[2], len = 1000)
    y.prior <- dPCPriPrec(x.prior, u = prior.u, a=prior.a, give_log = 0)
    
    tmb.post.median <- median(tmb.post.draws)
    inla.post.median <- median(inla.post.draws)
    param.name <- "clust prec"
    
    # ## test pc.prec.priors
    # prior.u <- c(.1, .25, .5, .75, 1)
    # prior.a <- c(.01, .05, .1, .25, .5)
    # par(mfrow=c(length(prior.u), length(prior.a)), mar= c(5, 4, 4, 2) + .1)
    # for(uu in prior.u){
    #   for(aa in prior.a){
    #     x <- seq(0, 200, len=1000)
    #     y <- dPCPriPrec(x, u=uu, a=aa)
    #     plot(x, y, main = sprintf('u = %0.3f || a = %0.3f', uu, aa))
    #     abline(v=100, col = 'red')
    #   }
    # }
    
    
  }

  ## get posterior samples (we'll use density curves)
  tmb.dens <- density(tmb.post.draws)
  inla.dens <- density(inla.post.draws)
  xrange <- xlim ## range(c(x.prior, tmb.dens$x, inla.dens$x))
  yrange <- range(c(y.prior, tmb.dens$y, inla.dens$y))

  prior.col <- "black"
  tmb.col   <- "red"
  inla.col  <- "blue"

  ## setup plot and plot prior
  plot(x.prior, y.prior, pch = ".", col = prior.col, main = param.name,
       xlim = xrange, ylim = yrange)
  lines(x.prior, y.prior, col = prior.col)

  ## plot tmb post and data
  lines(tmb.dens$x,tmb.dens$y, col = tmb.col)
  points(tmb.post.draws, rep(0, ndraws), col = alpha(tmb.col, 0.1), cex = 2, pch = '|')
  abline(v = tmb.post.median, col = tmb.col, lwd = 2)
    
  ## plot inla post and data
  lines(inla.dens$x,inla.dens$y, col = inla.col)
  points(inla.post.draws, rep(0, ndraws), col = alpha(inla.col, 0.1), cex = 2, pch = '|')
  abline(v = inla.post.median, col = inla.col, lwd = 2)

  ## plot truth
  abline(v = true.val, col = prior.col, lwd = 2)
  
  ## add legend
  legend("topright", legend = c("prior/truth", "tmb", "inla"), 
         col = c(prior.col, tmb.col, inla.col), lwd = rep(2, 3))
  
  ## and we also record the coverage for a parameter since we have the param draws here
  coverage_probs <- c(25, 50, 80, 90, 95)
  
  ## it's much faster to get all quantiles at once then to iteratively calculate them
  param.lui <- apply(rbind(inla.post.draws, tmb.post.draws), 1, quantile, 
               p = c((1 - coverage_probs/100)/2, 
                     coverage_probs/100 + (1 - coverage_probs/100) / 2), na.rm=T)
  
  ## add the binary coverage into the res table
  res[,type := c('truth', 'inla', 'tmb')]
  for(cc in 1:length(coverage_probs)){
    c <- coverage_probs[cc]
    li       <- param.lui[cc,]
    ui       <- param.lui[length(coverage_probs) + cc,]
    res[type=='inla', paste0(param, '.cov.',c) := (true.val >= li & true.val <= ui)[1]]
    res[type=='tmb', paste0(param, '.cov.',c) := (true.val >= li & true.val <= ui)[1]]
  }
  
}
dev.off()



## plot results in logit space or in prevalence space?
## TODO (?) plot.in.logit.space <- FALSE 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot comparisons of INLA and TMB fields
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting the comparison of INLA and TMB fields')

## randomly select pixels for plotting tmb v inla scatter
samp <- sample(cellIdx(ras_med_inla[[1]]),1e4) 

for(sum.meas in c('median','stdev')){ 
  
  if(sum.meas=='median'){
    
    if(data.lik == 'binom'){
      true  <- invlogit(true.rast[[1]])
      rinla <- ras_med_inla_p[[1]]
      rtmb  <- ras_med_tmb_p[[1]]
    }else if(data.lik == 'normal'){
      true <- true.rast[[1]]
      rinla <- ras_med_inla[[1]]
      rtmb  <- ras_med_tmb[[1]]
    }
    all.vec <- c(as.vector(rtmb), as.vector(rinla), as.vector(true))
    all.diff.vec <- as.vector(rinla - rtmb)
    diff.truth.zrange <- range(c( as.vector(rtmb-true), as.vector(rinla-true) ), na.rm=T)
    rast.list <- list('TRUE' = true,
                      'TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/experiment%06d_iter%06d_plot_03_median_rasters.png', out.dir, exp.lvid, exp.iter),
        height=12, width=12, units = 'in', res = 250)
  } ## sum.meas==median
  
  if(sum.meas=='stdev'){
    if(data.lik == 'binom'){
      rinla <- ras_sdv_inla_p[[1]]
      rtmb  <- ras_sdv_tmb_p[[1]]
    }else if(data.lik == 'normal'){
      rinla <- ras_sdv_inla[[1]]
      rtmb  <- ras_sdv_tmb[[1]]
    }
    all.vec <- c(as.vector(rtmb), as.vector(rinla))
    all.diff.vec <- as.vector(rinla - rtmb)
    rast.list <- list('TMB' = rtmb,
                      'INLA' = rinla)
    
    png(sprintf('%s/validation/experiment%06d_iter%06d_plot_04_stdev_rasters.png', out.dir, exp.lvid, exp.iter),
        height=8, width=8, units = 'in', res = 250)
  } ## sum.meas==stdev
  
  layout(matrix(1:length(rast.list) ^ 2, byrow = T, ncol = length(rast.list)))
  
  tmp <- subset(dt, period_id==1) ## for s-t

  ## get limits
  rast.zrange <- range(all.vec, na.rm = T)
  diff.zrange <- range(all.diff.vec, na.rm = T)
  
  ## set some plot args for mean and stdev separately
  if(sum.meas=='median'){
    ## set the legend.width to get nice plots
    lw <- 5
    tick.space.val  <- 0.5
    tick.space.truth.diff <- 0.25
    tick.space.diff <- 0.1
    mar <- c(0, 0, 1.4, 8)
  }
  if(sum.meas=='stdev'){
    ## set the legend.width to get nice plots
    lw <- 2
    tick.space.val  <- 0.1
    tick.space.diff <- 0.05
    mar <- c(0, 0, 1.4, 4)
  }
  
 

  for(i in 1:length(rast.list)){
    for(j in 1:length(rast.list)){

      if(i == j){
        ## plot raster
        par(mar = mar, bty='n')
        plot(rast.list[[i]],  maxpixel=1e7, col=rev(viridis(100)), axes=FALSE, legend.width = lw,
             axis.args=list(at=c(seq(rast.zrange[1], rast.zrange[2], by=tick.space.val), rast.zrange[2]),
                            labels=round(c(seq(rast.zrange[1], rast.zrange[2], by=tick.space.val), rast.zrange[2]), 3), 
                            cex.axis=0.6),
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1), main=paste0(names(rast.list)[i], ': ', sum.meas),
             zlim=rast.zrange)
      }

      if(i < j){
        ## plot scatter
        par(mar = c(4, 4, 2, 2),bty='n')
        plot(x=as.vector(rast.list[[j]])[samp],
             y=as.vector(rast.list[[i]])[samp],
             xlab='',
             ylab='',
             cex=.05, pch=19,
             main=paste0('(sub)SCATTER (', names(rast.list)[i], ' vs ', names(rast.list)[j], '): '))
        lines(x=rast.zrange, y=rast.zrange, col='red')
        title(xlab=names(rast.list)[j], line=-1, cex.lab=1.2)
        title(ylab=names(rast.list)[i], line=-1, cex.lab=1.2)
      }

      if(i > j){
        ## plot difference rasters
        
        if(sum.meas=='median' & j == 1){
          ## these are for truth - *_estiamte
          diff.col <- rev(plasma(100))
          dz  <- diff.truth.zrange
          tsd <- tick.space.truth.diff
          
        }else{
          ## these are for inla_est - tmb_est
          diff.col <- rev(cividis(100))
          dz  <- diff.zrange
          tsd <- tick.space.diff
        }
        par(mar = mar, bty='n')
        plot(rast.list[[i]] - rast.list[[j]],  maxpixel=1e7, col=diff.col, axes=FALSE, legend.width = lw,
             axis.args=list(at=sort(c(0, seq(dz[1], dz[2], by=tsd), dz[2])),
                            labels=round(sort(c(0, seq(dz[1], dz[2], by=tsd), dz[2])), 3),
                            cex.axis=0.6),
             legend.args=list(text='', side=2, font=1, line=0, cex=0.1),
             main=paste0('Diff: ', names(rast.list)[i], ' - ', names(rast.list)[j]),
             zlim=dz)
        points( x=tmp$long,y=tmp$lat, pch=19, cex=.1 )

      }
      
    } ## j
  } ## i

  dev.off()
  
} ## sum.meas

## ~~~~~~~~~~~~~~~~~~~~~~
## plot caterpillar plots
## ~~~~~~~~~~~~~~~~~~~~~~
message('------ plotting caterpillar plots')

png(sprintf('%s/validation/experiment%06d_iter%06d_plot_05_spatial_re_caterpillars.png', out.dir, exp.lvid, exp.iter),
      height=8, width=12, units = 'in', res = 250)

layout(matrix(1, 1, 1, byrow = TRUE))

## Compare mean and distribution of random effects
summ_gp_tmb  <- t(cbind((apply(epsilon_tmb_draws, 1, quantile,probs=c(.1,.5,.9)))))
summ_gp_inla <- t(cbind((apply(pred_s, 1, quantile,probs=c(.1,.5,.9)))))

## all time-space random effects
plot_d <- data.table(tmb_median = summ_gp_tmb[, 2],inla_median = summ_gp_inla[, 2],
                     tmb_low    = summ_gp_tmb[, 1],inla_low    = summ_gp_inla[, 1],
                     tmb_up     = summ_gp_tmb[, 3],inla_up     = summ_gp_inla[, 3])

plot_d$period <- factor(rep(1:nperiods, each=nrow(plot_d)/nperiods))
plot_d$loc    <- rep(1:(nrow(plot_d)/nperiods), rep=nperiods)
plot_d[, absdiff := abs(tmb_median-inla_median)]

## get the node locations and add them on
nodelocs <- mesh_s$loc
plot_d <- cbind(plot_d, nodelocs)

if(nrow(plot_d) > 2500)
  plot_d <- plot_d[sample(nrow(plot_d), 2500, replace=F), ]


## ## plot spatial RE differences
## ggplot(plot_d, aes(x=tmb_median,y=inla_median,col=period)) + theme_bw() +
## geom_point() + geom_line(aes(group=loc)) + geom_abline(intercept=0,slope=1,col='red') +
## ggtitle('Posterior Medians of Random Effects at Mesh Nodes, TMB v R-INLA. Connected dots same location different periods. ')

## plot locations where they are relatively different, are they near or far from data?
biggdiff <- plot_d[which(plot_d$absdiff > quantile(plot_d$absdiff, prob=0.80)), ]

plot(simple_polygon, main='Mesh nodes sized by abs difference TMB and R-INLA')
points(x=dt$long, y=dt$lat, pch=19, cex=(dt$N / max(dt$N))*.1)
points(x=plot_d$V1, y=plot_d$V2, pch=1, cex=plot_d$absdiff*5, col='red')

## catterpillar plot
plot_d <- plot_d[order(period,tmb_median)]
plot_d[ ,i := seq(1,.N), by = period]
gg_cat <- ggplot(plot_d, aes(i, tmb_median, col=i)) + theme_bw() + # [seq(1, nrow(plot_d), 5)]
          geom_linerange(aes(ymin = tmb_low, ymax = tmb_up), col='red', size=.8, alpha=.3) +
          geom_linerange(aes(x=i,ymin = inla_low, ymax = inla_up), col='blue', size=.8, alpha=.3) +
          facet_wrap(~period) +
          ggtitle('Comparison of random effects (10% to 90% quantiles) | BLUE == R-INLA | RED == TMB')
print(gg_cat)

dev.off()

## ###############################################
## 3) generate and summarize predictive metrics ##  
## ###############################################
message('------ making metrics to compare true and estimated surfaces')

## NOTE: these are in latent space for all models!
all.preds <- data.table(rbind(pred_tmb, pred_inla))

## make a data.table with prediction draws, model type, and truth
## NOTE! the truth and the median.fit are in logit-space!
non.na.idx <- which(!is.na(values(simple_raster)))
d <- data.table(truth        = rep(values(true.rast)[non.na.idx], 2),
                true.gp      = rep(values(true.gp)[non.na.idx], 2),
                model        = c(rep('tmb', nrow(pred_tmb)),
                                 rep('inla', nrow(pred_inla))), 
                median.fit   = apply(all.preds, 1, median))

## add on distance from nearest data observation to pixel

## first, get the centroids of the pixels
pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
pix.pts <- spTransform(pix.pts, CRS(geo.prj)) 
pix.pts@data <- data.frame(pix.pts@data, 
                           long = coordinates(pix.pts)[,1],
                           lat  = coordinates(pix.pts)[,2])
pix.pts.numeric <- as.data.frame(pix.pts@data)[,2:3]

## then calculate the distance from all the centroid pixels to 
## all the data observation locations, and take the smallest distance

min.obs.dist <- apply(raster::pointDistance(pix.pts.numeric, dt[,.(long, lat)],
                                      lonlat=T, allpairs = T),
                      1, min)/1000 ## distance comes out in meters

## break this into approx number of 5km pixels distant
min.obs.dist <- round(min.obs.dist / 5)
d[,obs.dist := rep(min.obs.dist, 2)]

## also, get the relative magnitude of the true gp
## break gp magnitude into deciles
d[,gp.dec := dplyr::ntile(true.gp, 10)]

# # ## make sure the decile is working
# par(mfrow=c(3, 4))
# plot(true.gp, main = 'true')
# for(i in 1:10){
#   plot(true.gp, main = sprintf('plotting decile %i', i))
#   ## message(sprintf('plotting decile %i', i))
#   dec.pts <- d[(.N/2+1):.N, which(gp.dec == i)]
#   ## points(pix.pts.numeric[dec.pts,], col=rev(brewer.pal(n=10, name='RdYlBu'))[i], pch=16, cex=.25)
#   points(pix.pts.numeric[dec.pts,], col=1, pch=16, cex=.25)
#   ## readline(prompt="Press [enter] to continue")
# }
# plot(true.gp, main='observations')
# points(dt[,.(long, lat)])

## get some coverage probs
coverage_probs <- c(25, 50, 80, 90, 95)
## it's much faster to get all quantiles at once then to iteratively calculate them
lui <- apply(all.preds, 1, quantile, 
             p = c((1 - coverage_probs/100)/2, 
                   coverage_probs/100 + (1 - coverage_probs/100) / 2), na.rm=T)
for(cc in 1:length(coverage_probs)){
  c <- coverage_probs[cc]
  message(paste0('-------- calcing ', c ,'% coverage of binomial prob'))
  li       <- lui[cc,]
  ui       <- lui[length(coverage_probs) + cc,]
  d[, paste0('pix.cov.',c)] <- d[['truth']] >= li & d[['truth']] <= ui
}

## get error and pred var
d[, error := truth - median.fit]
d[, var := apply(all.preds, 1, var)]

if(data.lik == 'binom'){
  ## NOTE: these are in latent space for all models!
  all.preds <- data.table(rbind(pred_tmb_p, pred_inla_p))

  ## make a data.table with prediction draws, model type, and truth
  ## NOTE! the truth and the median.fit are in logit-space!
  non.na.idx <- which(!is.na(values(simple_raster)))
  d[,`:=`(truth.p  = rep(values(true.rast.p)[non.na.idx], 2),
          model.p  = c(rep('tmb', nrow(pred_tmb_p)),
                       rep('inla', nrow(pred_inla_p))), 
          median.fit.p   = apply(all.preds, 1, median))]
  
  ## get some coverage probs
  coverage_probs <- c(25, 50, 80, 90, 95)
  ## it's much faster to get all quantiles at once then to iteratively calculate them
  lui <- apply(all.preds, 1, quantile, 
               p = c((1 - coverage_probs/100)/2, 
                     coverage_probs/100 + (1 - coverage_probs/100) / 2), na.rm=T)
  for(cc in 1:length(coverage_probs)){
    c <- coverage_probs[cc]
    message(paste0('-------- calcing ', c ,'% coverage of binomial prob'))
    li       <- lui[cc,]
    ui       <- lui[length(coverage_probs) + cc,]
    d[, paste0('p.pix.cov.',c)] = d[['truth.p']] >= li & d[['truth.p']] <= ui
  }

  ## get error and pred var
  d[, error.p := truth.p - median.fit.p]
  d[, var.p := apply(all.preds, 1, var)]
}

## MAKE COVERAGE PLOTS

require('tidyr')
message('------ making coverage plots')
cov.width <- 0.5 ## set the coverage interval width

## TODO add coverage plot of pixels from 7 
 
## TODO add coverage plot of params from 7

## make some plots of coverage by distance to observations
message('-------- making coverage vs dist to obs plots')

long.cov <- data.table(gather(d,
                              target_cov,  ## name of NEW key col
                              obs_cov,     ## name of NEW value col
                              pix.cov.25:pix.cov.95, ## cols to convert from wide to long
                              factor_key = TRUE))
## get the coverages calculated
nom.cov.names <- long.cov[,sort(unique(target_cov))] 
nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100
## and assign the numeric coverage for the rows
for(i in 1:length(nom.cov.names)){
  long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
}

## facet by observations
long.cov <- na.omit(long.cov)
long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                           groupvars = c('obs.dist', 'model', 'n_target_cov'),
                           conf.interval = cov.width)


## add on some relevant lvid params and save to use in the combination validation plot
long.cov.sum$lv.id     <- exp.lvid
long.cov.sum$n.clust   <- n.clust
long.cov.sum$nv        <- norm.var
long.cov.sum$clust.var <- clust.var
long.cov.sum$dl        <- data.lik
long.cov.sum$fit_type  <- long.cov.sum$model
## and save for validation
write.csv(long.cov.sum, sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv', out.dir, exp.lvid, exp.iter))

## make a plot of coverage as function of distance
## group the lines
long.cov.sum$line_group <- paste(long.cov.sum$model, long.cov.sum$obs.dist, sep='_')

pd <- position_dodge(0.05)
gg.dist.cov <- ggplot(long.cov.sum,
                      aes(x = n_target_cov, y = obs_cov,
                          shape = model, 
                          linetype = model,
                          color = factor(obs.dist),
                          group = line_group),
                      position=position_jitter(w=0.02, h=0.02)) + 
  #geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_abline(intercept = 0, slope = 1) +
  ##facet_wrap(. ~ n.clust, labeller = facet_labs) + 
  ggtitle('Coverage vs Distance') +
  ## fix the legends a bit
  labs(color = "Pixels to Nearest Obs", shape='Fit Type', linetype='Fit Type') + ## legend titles
  xlab('Nominal Coverage') + 
  ylab('Observed Monte Carlo Coverage')

if(interactive()){ ## then we can view
  print(gg.dist.cov)
}

ggsave(sprintf('%s/validation/experiment%06d_iter%06d_plot_06_coverage_distance_to_obs.png', out.dir, exp.lvid, exp.iter),
       plot = gg.dist.cov,
       device = 'png', units = 'in',
       width = 12, height = 8)

## make some plots of coverage by GP magnitude
message('-------- making coverage vs GP magnitude plots')

## facet by observations
long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
                           groupvars = c('gp.dec', 'model', 'n_target_cov'),
                           conf.interval = cov.width)

## add on some relevant lvid params and save to use in the combination validation plot
long.cov.sum$lv.id     <- exp.lvid
long.cov.sum$n.clust   <- n.clust
long.cov.sum$nv        <- norm.var
long.cov.sum$clust.var <- clust.var
long.cov.sum$dl        <- data.lik
long.cov.sum$fit_type  <- long.cov.sum$model
## and save for validation
write.csv(long.cov.sum, sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv', out.dir, exp.lvid, exp.iter))

## make a plot of coverage as function of distance
## group the lines
long.cov.sum$line_group <- paste(long.cov.sum$model, long.cov.sum$gp.dec, sep='_')

pd <- position_dodge(0.05)
gg.mag.cov <- ggplot(long.cov.sum,
                      aes(x = n_target_cov, y = obs_cov,
                          shape = model, 
                          linetype = model,
                          color = factor(gp.dec),
                          group = line_group),
                      position=position_jitter(w=0.02, h=0.02)) + 
  #geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_abline(intercept = 0, slope = 1) +
  ##facet_wrap(. ~ n.clust, labeller = facet_labs) + 
  ggtitle('Coverage vs True GP Magnitude') +
  ## fix the legends a bit
  labs(color = "GP Decile", shape='Fit Type', linetype='Fit Type') + ## legend titles
  xlab('Nominal Coverage') + 
  ylab('Observed Monte Carlo Coverage')

if(interactive()){ ## then we can view
  print(gg.mag.cov)
}

ggsave(sprintf('%s/validation/experiment%06d_iter%06d_plot_07_coverage_gpMagnitude.png', out.dir, exp.lvid, exp.iter),
       plot = gg.mag.cov,
       device = 'png', units = 'in',
       width = 12, height = 8)


message('-------- making final table of summary metrics and scores')
## summarize across pixels
surface.metrics <- data.table(cbind(mean.l      = d[, .(truth = mean(truth, na.rm = T)), by = c('model')], 
                                    mean.l.est  = d[, .(est = mean(median.fit, na.rm = T)), by = c('model')]$est, 
                                    bias        = d[, .(bias = mean(error, na.rm = T)), by = c('model')]$bias,
                                    rmse        = d[, .(rmse = sqrt(mean(error ^ 2, na.rm = T))), by = c('model')]$rmse,
                                    cor         = d[, .(cor = my.cor(truth, median.fit)), by = c('model')]$cor,
                                    CoV         = d[, .(cov = mean(error / var, na.rm = T)), by = c('model')]$cov,
                                    crps        = d[, .(crps = mean(crpsNormal(truth, median.fit, var), na.rm = T)), by = c('model')]$crps,
                                    pix.cov25       = d[, .(cov25 = mean(pix.cov.25, na.rm = T)), by = c('model')]$cov25,
                                    pix.cov50       = d[, .(cov50 = mean(pix.cov.50, na.rm = T)), by = c('model')]$cov50,
                                    pix.cov80       = d[, .(cov80 = mean(pix.cov.80, na.rm = T)), by = c('model')]$cov80,
                                    pix.cov90       = d[, .(cov90 = mean(pix.cov.90, na.rm = T)), by = c('model')]$cov90,
                                    pix.cov95       = d[, .(cov95 = mean(pix.cov.95, na.rm = T)), by = c('model')]$cov95))

if(data.lik == 'binom'){
  surface.metrics.p <- data.table(cbind(
                                    mean.p      = d[, .(truth.p = mean(truth.p, na.rm = T)), by = c('model')]$truth.p, 
                                    mean.p.est  = d[, .(est.p = mean(median.fit.p, na.rm = T)), by = c('model')]$est.p, 
                                    bias.p      = d[, .(bias.p = mean(error.p, na.rm = T)), by = c('model')]$bias.p,
                                    rmse.p      = d[, .(rmse.p = sqrt(mean(error.p ^ 2, na.rm = T))), by = c('model')]$rmse.p,
                                    cor.p       = d[, .(cor.p = my.cor(truth.p, median.fit.p)), by = c('model')]$cor.p,
                                    CoV.p       = d[, .(cov.p = mean(error.p / var.p, na.rm = T)), by = c('model')]$cov.p
                                    ))
  surface.metrics <- cbind(surface.metrics, surface.metrics.p)
}

## add together all two row summaries. first row inla, swcond row tmb
surface.metrics <- surface.metrics[order(mean.l.model), ]

## add on things from the res tab
res.addon <- t(rr[, c('R-INLA', 'TMB')])
colnames(res.addon) <- rr[, quantity]

## add on parameter coverage
pc.addon <- res[type != 'truth', .SD, .SDcols = names(res) %like% ".cov."]

## add on loopvars for easy summary later
lv.addon <- rbind(loopvars[exp.lvid, ], loopvars[exp.lvid, ])

## combine the pieces
summary.metrics <- cbind(surface.metrics, res.addon, pc.addon, lv.addon)

## note iteration
summary.metrics[, iter := exp.iter]

## save
write.csv(summary.metrics, sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv', out.dir, exp.lvid, exp.iter))
