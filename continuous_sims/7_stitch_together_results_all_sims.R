## NOTE: this script may be run in a clean env:

user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)
tmb_repo  <- "~/Documents/GitRepos/tmb_inla_comp"

## grab libraries and functions from MBG code
#setwd(core_repo)
#commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
#package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

## Load MBG packages and functions
#message('Loading in required R packages and MBG functions')
#source(paste0(core_repo, 'mbg_central/setup.R'))
#mbg_setup(package_list = package_list, repos = core_repo)

## load util funcs
## Now we can switch to the TMB repo
setwd(tmb_repo)
source('./realistic_sim_utils.R')

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(plyr)
library(ggplot2)
library(data.table)

## this script pulls together some overall results from an experiment run
# main.dir.name <- "2020_08_08_13_13_30" # using SPDE to sim GP
main.dir.name <- "spatial-stats-sub-2020_09_05_21_56_12" # using RandomFields to sim GP
main.dir.root <- "/ihme/scratch/users/azimmer/tmb_inla_sim"
main.dir.root <- "/home/merms/Documents/Research/2020_tmb_v_inla/tmb_inla_sim"

main.dir <- file.path(main.dir.root, main.dir.name)
compar.dir <- file.path(main.dir, "comparisons")
# compar.dir <- "/ihme/homes/azimmer/tmb_inla_sim/2020_09_05_21_56_12/comparisons"
# dir.create(compar.dir, recursive = T, showWarnings = F)

## read in all experiment parameters
loopvars <- fread(file = paste0(main.dir, '/loopvars.csv'), stringsAsFactors = F)

## print columns of loopvar that vary
## so we can easily see what's going on in the experiments that fail...
loopvars[, !apply(loopvars, MARGIN = 2,
                  FUN=function(x){col.var <- sort(x, decreasing=F)[1] == sort(x, decreasing=T)[1]
                    if(is.na(col.var)){
                      return(TRUE)
                    }else{
                      return(col.var)}
                  }), with=F]

## if we haven't already done this,
## read in the summary metrics log from each iteration of each experiment
missing.files <- data.table(lvid=integer(),
                            iter=integer())
if(!file.exists(sprintf('%sall_summary_metrics.csv', compar.dir))){

  summary.metrics.list <- cov.gp.list <- cov.dist.list <- list(NULL)
  list.ind <- 1

  for(lvid in 1:nrow(loopvars)){
    message(sprintf('--loading in summary metrics from %i of %i', lvid, nrow(loopvars)))

    out.dir  <- sprintf('%s/%06d', main.dir, lvid)

    for(iter in 1:loopvars$n.sim[1]){ ## n.sim is the same for all experiments

      ## check to see that the file is there
      if(!file.exists(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv',
                              out.dir, lvid, iter))){
        message(sprintf('----WARNING!! summary file for exp %06d, iter %06d does not exist!', lvid, iter))
        missing.files <- rbind(missing.files,
                               data.table(lvid=lvid,
                                          iter=iter))
      }else{

        summary.metrics.list[[list.ind]] <- fread(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv',
                                                          out.dir, lvid, iter))[,lvid:=lvid]

        cov.gp.list[[list.ind]]          <- fread(sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv',
                                                          out.dir, lvid, iter))[,lvid:=lvid][,iter:=iter]

        cov.dist.list[[list.ind]]        <- fread(sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv',
                                                          out.dir, lvid, iter))[,lvid:=lvid][,iter:=iter]

        list.ind <- list.ind +  1

        # if(lvid==1 & iter==1){
        #   summary.metrics <- fread(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv',
        #                                    out.dir, lvid, iter))[,lvid:=lvid]
        #   cov.gp   <- fread(sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv',
        #                             out.dir, lvid, iter))[,lvid:=lvid]
        #   cov.dist <- fread(sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv',
        #                             out.dir, lvid, iter))[,lvid:=lvid]
        # }else{
        #   summary.metrics <- rbind(summary.metrics,
        #                            fread(sprintf('%s/validation/experiment%06d_iter%06d_summary_metrics.csv',
        #                                          out.dir, lvid, iter))[,lvid:=lvid], fill=T)
        #   cov.gp   <- rbind(cov.gp,
        #                     fread(sprintf('%s/validation/experiment%06d_iter%06d_GP_magnitude_coverage_summary.csv',
        #                                   out.dir, lvid, iter))[,lvid:=lvid], fill = T)
        #   cov.dist <- rbind(cov.dist,
        #                     fread(sprintf('%s/validation/experiment%06d_iter%06d_distance_coverage_summary.csv',
        #                                   out.dir, lvid, iter))[,lvid:=lvid], fill = T)
        # }

      } ## file.exists check
    }   ## loading all iterations within
  }     ## all experiments (ie loopvars loop)

  ## convert all columns of summary metrics to character to ensure constant typing of each column,
  ##  then rbind, then convert columns to appropriate types (numerics, logicals, etc)
  summary.metrics.list2 <- lapply(summary.metrics.list,
                                  function(x){
                                    # make sure each column only shows up once
                                    x <- x[, .SD, .SDcols = unique(names(x))]
                                    # convert to characters
                                    x[, names(x) := lapply(.SD, as.character)]
                                  })
  ## combine the lists
  summary.metrics <- rbindlist(summary.metrics.list2, fill=T)
  cov.gp          <- rbindlist(cov.gp.list, fill=T)
  cov.dist        <- rbindlist(cov.dist.list, fill=T)

  ## save the combined metrics for the study
  write.csv(summary.metrics, file = sprintf('%sall_summary_metrics.csv', compar.dir))
  write.csv(cov.gp, file = sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
  write.csv(cov.dist, file = sprintf('%sall_dist_coverage_metrics.csv', compar.dir))

  if(nrow(missing.files) > 0){
    message('-- WARNING! not all experiments and iters had summary metrics files. missing:')
    print(missing.files)
  }else{
    message('--summary metrics from all experiments and all iterations succesfully combined and saved')
  }

}else{

  ## reload the prepped file
  message('--reloading in combined summary metrics')
  summary.metrics <- fread(sprintf('%sall_summary_metrics.csv', compar.dir))
}

message('summary metrics from all experiments and all iterations prepped and ready to go')

## convert columns of summary.metrics to appropriate type
numCols <- c( "mean.l.truth",
             "mean.l.est"           ,"bias",
             "rmse"                 ,"cor",
             "CoV"                  ,"crps",
             "pix.cov25"            ,"pix.cov50",
             "pix.cov80"            ,"pix.cov90",
             "pix.cov95"            ,
             "st_mesh_nodes",
             "cores"                ,"s_mesh_max_edge",
             "s_mesh_cutoff"        ,"draws",
             "fit_time"             ,"pred_time",
             "pt_tmb_sdreport_time" ,"pt_get_draws_time",
             "convergence.fails",
             "fe_int_mean"          ,"fe_int_sd",
             "gauss_var_mean"       ,"gauss_var_sd",
             "matern_range_mean"    ,"matern_range_sd",
             "matern_sigma_mean"    ,"matern_sigma_sd",
             "gauss_prec_mean"      ,
             "year_list"            ,
             "alpha"                ,"sp.range",
             "sp.var"               ,"sp.alpha",
             "clust.var"            ,"t.rho",
             "n.clust",
             "m.clust"              ,
             "ndraws"               ,"n.sim",
             "norm.var"             ,"iter",
             "lvid"                 ,"fe_access2_mean",
             "fe_mapincidence_mean" ,"fe_access2_sd",
             "fe_mapincidence_sd"   ,
             "clust_var_mean"       ,"clust_var_sd",
             "clust_prec_mean"      ,
             "mean.p"              ,"mean.p.est",
             "bias.p"              ,"rmse.p",
             "cor.p"               ,"CoV.p")

# logical columns
logCols <- c(
  "convergence"          ,"alpha.cov.25",
  "alpha.cov.50"         ,"alpha.cov.80",
  "alpha.cov.90"         ,"alpha.cov.95",
  "gauss.prec.cov.25"    ,"gauss.prec.cov.50",
  "gauss.prec.cov.80"    ,"gauss.prec.cov.90",
  "gauss.prec.cov.95"    ,"sp.range.cov.25",
  "sp.range.cov.50"      ,"sp.range.cov.80",
  "sp.range.cov.90"      ,"sp.range.cov.95",
  "sp.sigma.cov.25"      ,"sp.sigma.cov.50",
  "sp.sigma.cov.80"      ,"sp.sigma.cov.90",
  "sp.sigma.cov.95"      ,"bias.correct",
  "sd.correct"           ,
  "fix.locs"             ,"fix.gp",
  "beta.cov.25",
  "beta.cov.50"          ,"beta.cov.80",
  "beta.cov.90"          ,"beta.cov.95",
  "clust.prec.cov.25",
  "clust.prec.cov.50"   ,"clust.prec.cov.80",
  "clust.prec.cov.90"   ,"clust.prec.cov.95")

# convert
summary.metrics[,(numCols):= lapply(.SD, as.numeric), .SDcols = numCols]
summary.metrics[,(logCols):= lapply(.SD, as.logical), .SDcols = logCols]

## process a few things and make some necessary columns for plotting

## set the order of clust.var so that NA -> None and comes first
## set factor order of clust.var
clust.var <- as.character(summary.metrics$clust.var)
non.na.cv <- as.character(sort(unique(na.omit(summary.metrics$clust.var))))
clust.var[is.na(clust.var)] <- 'None'
summary.metrics$clust.var.cat <- factor(clust.var, levels = c('None', non.na.cv))

## get the true logkappa and logtau params
summary.metrics[,logkappa := log(sqrt(8) / sp.range)]
summary.metrics[,logtau   := log(sqrt(1 / (4 * pi * exp(logkappa) ^ 2 * sp.var)))]

## make new labels for INLA_EB, INLA_CCD, TMB
summary.metrics[inla.int.strat == 'eb' & mean.l.model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_EB_G']
summary.metrics[inla.int.strat == 'eb' & mean.l.model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_EB_S']
summary.metrics[inla.int.strat == 'eb' & mean.l.model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_EB_L']
summary.metrics[inla.int.strat == 'ccd' & mean.l.model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_CCD_G']
summary.metrics[inla.int.strat == 'ccd' & mean.l.model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_CCD_S']
summary.metrics[inla.int.strat == 'ccd' & mean.l.model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_CCD_L']
summary.metrics[mean.l.model == 'tmb', fit_type := 'TMB']

## ## and do this for cov.gp
## moved to where we need these objects for memory reasons
## cov.gp          <- fread(sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
## to.merge <- unique(summary.metrics[,.(lvid, inla.int.strat, inla.approx, st_mesh_nodes)])
## cov.gp[,V1:=NULL]
## cov.gp <- merge(cov.gp, to.merge, by='lvid', sort=FALSE)

## cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_EB_G']
## cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_EB_S']
## cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_EB_L']
## cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_CCD_G']
## cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_CCD_S']
## cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_CCD_L']
## cov.gp[model == 'tmb', fit_type := 'TMB']

## ## and for cov.dist
## cov.dist        <- fread(sprintf('%sall_dist_coverage_metrics.csv', compar.dir))
## to.merge <- unique(summary.metrics[,.(lvid, inla.int.strat, inla.approx, st_mesh_nodes)])
## cov.dist[,V1:=NULL]
## cov.dist <- merge(cov.dist, to.merge, by='lvid', sort=FALSE)

## cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_EB_G']
## cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_EB_S']
## cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_EB_L']
## cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_CCD_G']
## cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_CCD_S']
## cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_CCD_L']
## cov.dist[model == 'tmb', fit_type := 'TMB']

## ####################################################################
## ####################################################################
## now that should be everything!
## we can make a bunch of plots
## after we format the data into long for ggplot
## ####################################################################
## ####################################################################

## Gather columns into key-value pairs and move them from wide to long format
## for plotting pixel level covereage
long.cov <- data.table(gather(summary.metrics,
                              target_cov,  ## name of NEW key col
                              obs_cov,     ## name of NEW value col
                              pix.cov25:pix.cov95, ## cols to convert from wide to long
                              factor_key = TRUE
                              ))

## get the coverages calculated
nom.cov.names <- long.cov[,sort(unique(target_cov))]
nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100

## and assign the numeric coverage for the rows
for(i in 1:length(nom.cov.names)){
  long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
}

## calculate noise to spatial ratio
long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0

## process the total fitting and predict times
long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
long.cov[, total.time :=  pred.time + fit.time]

## set cluster var NA to 0. NOTE: this is fit w/o cluster var params!
long.cov[is.na(clust.var), clust.var:=0]

## ##########################################################################
## ##########################################################################
## PIXEL COVERAGE: facets by number of observations
## ##########################################################################
## ##########################################################################

## loop through binomial and normal with different obs variances
## make one plot for each different data model

message('plotting average pixel coverage plots by sample size')

loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                          nv=c(NA, loopvars[, unique(norm.var)]))

cov.widths <- c(0.5, 0.8, 1.0) ## set the coverage interval width

cov.widths <- 0.8
compar.dir <- "/home/merms/Dropbox/AaronJon/TMB_SpatialStats/figures/cont_sim_res/"

for(cov.width in cov.widths){
  for(dl in c('binom', 'normal')){
    for(inla.int in c('EB_G', "EB_S", "EB_L", 'CCD_G', "CCD_S", "CCD_L")){
      for(with.tmb in "both"){ #}c("only", "not", "both")){ # plot TMB only, INLA only, both [on the same plots]

        # we only need to plot the TMB only one once, not once for every INLA!
        if(with.tmb=="only" & inla.int != "EB_G"){
          next
        }

        message(sprintf('-- plotting for dl=%s and inla.int=%s and with.tmb=%s with %0.0f%% ints',
                        dl, inla.int, with.tmb, cov.width*100))

        ## make the plot title
        pt <- sprintf('Average Pixel Coverage across Spatial Domain for %s Observations:\n',
                      stringr::str_to_title(dl))
        pt <- paste(pt, sprintf('Median and %0.0f%% Quantile Bars across Monte Carlo Replicates', cov.width*100))

        # if(dl == 'binom') {
        #   long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
        # }else if(dl == 'normal') {
        #   long.cov.sub <- subset(long.cov, data.lik == dl)
        # }

        # subset to particular liklihood and inla hyper param integration type
        if(with.tmb=="only"){
          ## only plot TMB
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   grepl("TMB", long.cov$fit_type))
        }
        if(with.tmb=="not"){
          # only plot INLA results
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   grepl(inla.int, long.cov$fit_type))

        }
        if(with.tmb=="both"){
          # plot INLA and TMB results
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   (grepl(inla.int, long.cov$fit_type)  |
                                      grepl("TMB", long.cov$fit_type)))
        }

        long.cov.sub[fit_type == 'INLA_CCD_L', fit_type := 'INLA: CCD + L']
        long.cov.sub[fit_type == 'INLA_EB_S', fit_type := 'INLA: EB + S']

        ## facet by observations
        long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                                   groupvars = c('n_target_cov', 'fit_type', 'n.clust', 'data.lik', 'clust.var.cat', 'norm.var'),
                                   conf.interval = cov.width)
        ## groups to plot lines
        long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$clust.var, sep='_')

        ## set facet names using the labeller
        ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
        ## facet_labs <- as_labeller(eval(parse(text=paste0('c(', paste0('`', sort(unique(long.cov$n.clust)), '`', '=', '\'', paste0('Num. Cluster Obs: ', sort(unique(long.cov$n.clust))), '\'', sep='', collapse=','), ')'))))

        ## rename for plotting
        setnames(long.cov.sum, 'n.clust',  'Number of Clusters')
        setnames(long.cov.sum, 'norm.var', 'Observation Variance')

        pd <- position_dodge(0.05)
        fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                          aes(x = n_target_cov, y = obs_cov,
                                              shape = fit_type,
                                              linetype = fit_type,
                                              color = clust.var.cat,
                                              group = line_group),
                                          position=position_jitter(w=0.02, h=0.02)) +
          geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
          geom_line(position=pd) +
          geom_point(position=pd) +
          geom_abline(intercept = 0, slope = 1) +
          ggtitle(pt) +
          ## fix the legends a bit
          labs(color = "Cluster Var.", shape='Method', linetype='Method') + ## legend titles
          theme(panel.spacing.x = unit(4, "mm")) +
          xlab(expression(Target~Coverage:~alpha)) +
          ylab('Mean Coverage of Field Estimates')

        if(dl=='binom'){
          fit_coverage_CI_summary <- fit_coverage_CI_summary +
            facet_wrap(. ~ `Number of Clusters`, labeller = label_both)
          plot.w = 10
          plot.h = 7
        }
        if(dl == 'normal'){
          fit_coverage_CI_summary <- fit_coverage_CI_summary +
            facet_grid(`Number of Clusters` ~ `Observation Variance`, labeller = label_both)
          plot.w = 10
          plot.h = 12
        }

        ## if(interactive()){ ## then we can view
        ##   print(fit_coverage_CI_summary)
        ## }

        ggsave(sprintf('%s/plot_01_%s_pixel_coverage_summary_by_observations_%s%s%s_%i_intervals.png',
                       compar.dir,
                       dl,
                       ifelse(with.tmb %in% c("only", "both"), "_TMB", ""),
                       ifelse(with.tmb %in% c("both"), "_and", ""),
                       ifelse(with.tmb != "only", sprintf("_INLA_%s", inla.int), ""),
                       cov.width*100),
               plot = fit_coverage_CI_summary,
               device = 'png', units = 'in',
               width = plot.w, height = plot.h)
      }
    }
  }
}


## ##########################################################################
## ##########################################################################
## PIXEL COVERAGE: facets by noise to spatial signal
## ##########################################################################
## ##########################################################################

message('plotting average pixel coverage plots by spatial noise:signal ratio')

loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
                          nv=c(NA, loopvars[, unique(norm.var)]))

cov.widths <- c(0.5, 0.8, 1.0) ## set the coverage interval width
for(cov.width in cov.widths){
  for(dl in c('binom', 'normal')){
    for(inla.int in c('EB_G', "EB_S", "EB_L", 'CCD_G', "CCD_S", "CCD_L")){
      for(with.tmb in c("only", "not", "both")){ # plot TMB only, INLA only, both [on the same plots]

        # we only need to plot the TMB only one once, not once for every INLA!
        if(with.tmb=="only" & inla.int != "EB_G"){
          next
        }

        message(sprintf('-- plotting for dl=%s and inla.int=%s and with.tmb=%s with %0.0f%% ints',
                        dl, inla.int, with.tmb, cov.width*100))

        ## make the plot title
        pt <- sprintf('Average Pixel Coverage given %s Observations',
                      stringr::str_to_title(dl))
        pt <- paste(pt, sprintf('\n %0.0f%% Intervals of Estimates across Simulations', cov.width*100))

        # if(dl == 'binom') {
        #   long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
        # }else if(dl == 'normal') {
        #   long.cov.sub <- subset(long.cov, data.lik == dl)
        # }

        # subset to particular liklihood and inla hyper param integration type
        if(with.tmb=="only"){
          ## only plot TMB
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   grepl("TMB", long.cov$fit_type))
        }
        if(with.tmb=="not"){
          # only plot INLA results
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   grepl(inla.int, long.cov$fit_type))

        }
        if(with.tmb=="both"){
          # plot INLA and TMB results
          long.cov.sub <- subset(long.cov,
                                 data.lik == dl &
                                   (grepl(inla.int, long.cov$fit_type)  |
                                      grepl("TMB", long.cov$fit_type)))
        }

        ## facet by noise:spatial var ratio
        long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                                   groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
                                               'noise_spatial_ratio', 'norm.var'),
                                   conf.interval=cov.width)

        ## groups to plot lines
        long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')

        ## rename for plotting
        setnames(long.cov.sum, 'noise_spatial_ratio', 'Cluster Var over Sp Var')
        setnames(long.cov.sum, 'norm.var', 'Obs Var')


        # ## set facet names using the labeller
        # ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
        # clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)]
        # ## sort by unique noise_spatial_ratio order
        # clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
        # facet_labs <- as_labeller(eval(parse(text =
        #                                        paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'',
        #                                                            sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f',
        #                                                                    clust.sp.pairs$clust.var,
        #                                                                    clust.sp.pairs$sp.var,
        #                                                                    sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
        #                                               ')'))))

        pd <- position_dodge(0.05)
        fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                          aes(x=n_target_cov, y=obs_cov,
                                              shape = fit_type,
                                              linetype = fit_type,
                                              colour=factor(n.clust),
                                              group = line_group),
                                          position=position_jitter(w=0.01, h=0.01)) +
          geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
          geom_line(position=pd) +
          geom_point(position=pd) +
          geom_abline(intercept = 0, slope = 1) +
          ## facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
          ggtitle(pt) +
          ## fix the legends a bit
          labs(color = 'Num. Obs', shape='Method', linetype = 'Method') + ## legend titles
          xlab(expression(Target~Coverage:~alpha)) +
          ylab('Mean Coverage of Field Estimates')

        if(dl=='binom'){
          fit_coverage_CI_summary <- fit_coverage_CI_summary +
            facet_wrap(. ~ `Cluster Var over Sp Var`, labeller = label_both)
        }
        if(dl == 'normal'){
          fit_coverage_CI_summary <- fit_coverage_CI_summary +
            facet_grid(`Cluster Var over Sp Var` ~ `Obs Var`, labeller = label_both)
        }

        if(interactive()){ ## then we can view
          print(fit_coverage_CI_summary)
        }

        ggsave(sprintf('%s/plot_02_%s_pixel_coverage_summary_by_noiseSpatialRatio_%s%s%s_%i_intervals.png',
                       compar.dir,
                       dl,
                       ifelse(with.tmb %in% c("only", "both"), "_TMB", ""),
                       ifelse(with.tmb %in% c("both"), "_and", ""),
                       ifelse(with.tmb != "only", sprintf("_INLA_%s", inla.int), ""),
                       cov.width*100),
               plot = fit_coverage_CI_summary,
               device = 'png', units = 'in',
               width = 12, height = 8)
      }
    }
  }
}


## ##########################################################################
## ##########################################################################
## effects bias & coverage
## ##########################################################################
## ##########################################################################

message('plotting effect bias')

## get all the fixed effects
non.pix.eff.sd.colnames <- names(summary.metrics)[grep('_sd$', names(summary.metrics))]
## now, get the first part of the name
non.pix.eff.b.colnames <- paste0(unlist(lapply( non.pix.eff.sd.colnames,
                                               function(x){
                                                 substr(x, 1, nchar(x)-3)
                                               })), '_bias')

(summary.metrics[,`Intercept` := fe_int_mean - alpha])
(summary.metrics[,`Obs. Var.` := gauss_var_mean - norm.var])
(summary.metrics[,`Matern Range` := matern_range_mean - sp.range])
(summary.metrics[,`Matern SD` := matern_sigma_mean - sqrt(sp.var)])
(summary.metrics[,`Cluster Var.` := clust_var_mean - clust.var])
(summary.metrics[,`Access` := fe_access2_mean - eval(parse(text=betas))[1]])
# NOTE! assumes access2 is first and mapincidence is second! count be written better....
(summary.metrics[,`Malaria` := fe_mapincidence_mean - eval(parse(text=betas))[2]])

non.pix.eff.b.colnames <- c('Intercept',
                            'Access',
                            'Malaria',
                            'Obs. Var.',
                            'Cluster Var.',
                            'Matern SD',
                            'Matern Range')

## get some coverage probs
key.summ.metrics <- summary.metrics[,c(non.pix.eff.b.colnames,
                                       'data.lik', 'norm.var', 'n.clust',
                                       'fit_type', 'clust.var.cat',
                                       "inla.int.strat", "inla.approx"), with=F]
fe.mean.long <- melt(key.summ.metrics, id=c('data.lik', 'norm.var', 'n.clust',
                                            'fit_type', 'clust.var.cat',
                                            "inla.int.strat", "inla.approx"))

## drop NAs. eg there is no gauss_prec_bias in a binom model,
## and there is no clust_prec_bias when clust.var==NA
fe.mean.long <- fe.mean.long[!is.na(value),]

cov.widths <- c(0.5, 0.8, 1.0) ## set the coverage interval width
cov.widths <- 0.8
compar.dir <- "/home/merms/Dropbox/AaronJon/TMB_SpatialStats/figures/cont_sim_res/"

for(cov.width in cov.widths){
  for(dl in c('binom', 'normal')){
    for(inla.int in c('EB_G', "EB_S", "EB_L", 'CCD_G', "CCD_S", "CCD_L")){
      for(with.tmb in "both"){ #}c("only", "not", "both")){ # plot TMB only, INLA only, both [on the same plots]

        # we only need to plot the TMB only one once, not once for every INLA!
        if(with.tmb=="only" & inla.int != "EB_G"){
          next
        }

        message(sprintf('-- plotting for dl=%s and inla.int=%s and with.tmb=%s with %0.0f%% ints',
                        dl, inla.int, with.tmb, cov.width*100))

        ## make the plot title
        pt <- sprintf('Parameter Bias for %s Observations:',
                      stringr::str_to_title(dl))
        pt <- paste(pt, sprintf('Median and %0.0f%% Quantile Bars across Monte Carlo Replicates', cov.width*100))

        # if(dl == 'binom') {
        #   fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
        # }else if(dl == 'normal') {
        #   fe.mean.long.sub <- subset(fe.mean.long, data.lik == dl)
        # }

        # subset to particular liklihood and inla hyper param integration type
        if(with.tmb=="only"){
          ## only plot TMB
          long.cov.sub <- subset(fe.mean.long,
                                 data.lik == dl &
                                   grepl("TMB", fe.mean.long$fit_type))
        }
        if(with.tmb=="not"){
          # only plot INLA results
          long.cov.sub <- subset(fe.mean.long,
                                 data.lik == dl &
                                   grepl(inla.int, fe.mean.long$fit_type))

        }
        if(with.tmb=="both"){
          # plot INLA and TMB results
          long.cov.sub <- subset(fe.mean.long,
                                 data.lik == dl &
                                   (grepl(inla.int, fe.mean.long$fit_type)  |
                                      grepl("TMB", fe.mean.long$fit_type)))
        }

        long.cov.sub[fit_type == 'INLA_CCD_L', fit_type := 'INLA: CCD + L']
        long.cov.sub[fit_type == 'INLA_EB_S', fit_type := 'INLA: EB + S']

        ## facet by ???
        fe.mean.long.sum <-  summarySE(long.cov.sub, measurevar="value",
                                       groupvars=c("data.lik", 'norm.var', 'n.clust',
                                                   "fit_type", 'clust.var.cat',
                                                   'variable'),
                                       conf.interval = cov.width)

        ## groups to plot lines
        ## long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')

        ## rename for plotting
        setnames(fe.mean.long.sum, 'n.clust',  'Number of Clusters')
        setnames(fe.mean.long.sum, 'norm.var', 'Observation Variance')

        pd <- position_dodge(.25)
        fe_bias_summary <- ggplot(fe.mean.long.sum,
                                  aes(x=log(`Number of Clusters`), y=med,
                                      shape = fit_type,
                                      linetype = fit_type,
                                      colour = clust.var.cat)) +
          geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
          geom_line(position=pd) +
          geom_point(position=pd, size=2) +
          geom_abline(intercept = 0, slope=0) +
          ## facet_wrap(. ~ variable, scales='free_y') +
          ggtitle(pt) +
          ## fix the legends a bit
          labs(color = 'Cluster Var.', shape='Method', linetype = 'Method') + ## legend titles
          xlab('Number of Clusters') +
          scale_x_continuous(breaks = log(sort(fe.mean.long.sum$`Number of Clusters`)), ## add log x-axis labels
                             labels = paste0('ln(', sort(fe.mean.long.sum$`Number of Clusters`), ')')) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', family='sans')) + ## rotate x-axis ticks
          ylab('Bias: (Estimate - True)')

        if(dl=='binom'){
          fe_bias_summary <- fe_bias_summary +
            facet_wrap(. ~ variable,
                       labeller = labeller(.cols = label_value),
                       scales='free_y')
          plot.w <- 10
          plot.h <- 7
        }
        if(dl == 'normal'){
          fe_bias_summary <- fe_bias_summary +
            facet_grid(variable ~ `Observation Variance`,
                       labeller = labeller(.rows = label_value, .cols = label_both),
                       scales='free')
          plot.w <- 10
          plot.h <- 12
        }

        ## if(interactive()){ ## then we can view
        ##   print(fe_bias_summary)
        ## }

        ggsave(sprintf('%s/plot_03_%s_param_bias_summary_%s%s%s_%i_intervals.png',
                       compar.dir,
                       dl,
                       ifelse(with.tmb %in% c("only", "both"), "_TMB", ""),
                       ifelse(with.tmb %in% c("both"), "_and", ""),
                       ifelse(with.tmb != "only", sprintf("_INLA_%s", inla.int), ""),
                       cov.width*100),
               plot = fe_bias_summary,
               device = 'png', units = 'in',
               width = plot.w, height = plot.h)
      }
    }
  }
}


## #####################################
## #####################################
## TIMING PLOTS
## #####################################
## #####################################

# 25 distinct colors
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)

summary.metrics[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
summary.metrics[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
summary.metrics[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
summary.metrics[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
summary.metrics[, total.time :=  pred.time + fit.time]

long.cov <- data.table(gather(summary.metrics,
                              operation,
                              time_s,
                              fit.time:total.time,
                              factor_key = TRUE))

# change names for plotting
long.cov[operation == "fit.time", operation := "Fit"]
long.cov[operation == "pred.time", operation := "Prediction"]
long.cov[operation == "total.time", operation := "Total"]

## facet by observations
long.cov.sub = copy(long.cov)
long.cov.sub[fit_type == 'INLA_CCD_L', fit_type := 'INLA: CCD + L']
long.cov.sub[fit_type == 'INLA_CCD_G', fit_type := 'INLA: CCD + G']
long.cov.sub[fit_type == 'INLA_CCD_S', fit_type := 'INLA: CCD + S']
long.cov.sub[fit_type == 'INLA_EB_S', fit_type := 'INLA: EB + S']
long.cov.sub[fit_type == 'INLA_EB_L', fit_type := 'INLA: EB + L']
long.cov.sub[fit_type == 'INLA_EB_G', fit_type := 'INLA: EB + G']
long.cov.sub[data.lik == 'binom', data.lik := 'Binom']
long.cov.sub[data.lik == 'normal', data.lik := 'Normal']

long.cov.sum <-  summarySE(long.cov.sub, measurevar="time_s",
                           groupvars=c("operation","fit_type", 'n.clust', 'st_mesh_nodes', 'data.lik'), conf.interval = 0.8)
setnames(long.cov.sum, 'st_mesh_nodes', 'Number of Spatial REs')
setnames(long.cov.sum, 'operation', 'Operation')

long.cov.sum$line_group <- paste(  long.cov.sum$fit_type,   long.cov.sum$data.lik)

pd <- position_dodge(.25)
fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                  aes(x=log(n.clust), y=time_s, colour=fit_type, group = line_group, linetype = data.lik)) +
  geom_errorbar(aes(ymin=l.ci, ymax=u.ci), width=.01, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  scale_x_continuous(breaks = log(sort(long.cov.sum$n.clust)), ## add log x-axis labels
                     labels = paste0('ln(', sort(long.cov.sum$n.clust), ')')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', family='sans')) + ## rotate x-axis ticks
  facet_grid(`Number of Spatial REs` ~ Operation,
             labeller = labeller(.rows = label_both, .cols = label_value),
             scales='free_y') +
  ggtitle('Comparison of Fit and Predict Times: Median and 80% Quantile Bars across Monte Carlo Replicates') +
  ylab("Time (seconds)") +
  xlab("Number of Clusters") +
  labs(colour="Method", linetype = 'Obs. Type') +
  scale_colour_manual(values = c25[c(20, 2, 3, 4, 5, 7, 6)])

if(interactive()){
  ## then we can view
  print(fit_coverage_CI_summary)
}

ggsave(sprintf('%s/fit_time_num_res.png', compar.dir),
       plot = fit_coverage_CI_summary,
       device = 'png', units = 'in',
       width = 9, height = 9)




## fit_coverage_CI_summary <- ggplot(long.cov.sum,
##                                   aes(x=n.clust, y=time_s, colour=fit_type, group = fit_type)) +
##   geom_errorbar(aes(ymin=l.ci, ymax=u.ci), width=.01) +
##   geom_line() +
##   geom_point() +
##   facet_wrap(operation ~ n.clust) + ggtitle(sprintf('Comparison of fit time (sec) in: %s', loopvars$data.lik[1]))

## if(interactive()){
##   ## then we can view
##   print(fit_coverage_CI_summary)
## }

## ggsave(sprintf('%s/%s_fit_time_nclust.png', compar.dir, loopvars$data.lik[1]),
##        plot = fit_coverage_CI_summary,
##        device = 'png', units = 'in',
##        width = 12, height = 12)


## ##########################################################################
## ##########################################################################
## PLOT COVERAGE BY GP DECILE
## ##########################################################################
## ##########################################################################

message('plotting coverage by gp decile')

## load releavant data and merge on gp decile coverage
to.merge <- unique(summary.metrics[,.(lvid, inla.int.strat, inla.approx, st_mesh_nodes)])

## clean up big mem objects
rm(long.cov); rm(summary.metrics); rm(fe.mean.long); rm(long.cov.sub); rm(key.summ.metrics); rm(clust.var)
gc()

compar.dir <- file.path(main.dir, "comparisons")
cov.gp          <- fread(sprintf('%sall_gp_coverage_metrics.csv', compar.dir))
cov.gp[,V1:=NULL]
cov.gp <- merge(cov.gp, to.merge, by='lvid', sort=FALSE)
for(i in 1:3){gc()}

cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_EB_G']
cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_EB_S']
cov.gp[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_EB_L']
cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_CCD_G']
cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_CCD_S']
cov.gp[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_CCD_L']
cov.gp[model == 'tmb', fit_type := 'TMB']

## cov.gp is already long, and mostly ready for plotting
## we need to average over the right things and make a column for line groupings

loop.params <- data.table(expand.grid(dl=c('binom', 'normal'),
                                      nv=c(NA, loopvars[, unique(norm.var)]),
                                      inla.strat = c("CCD_S"), #c('EB_G', "EB_S",  "EB_L",
                                      #  'CCD_G', "CCD_S",
                                      #                                                     "CCD_L"),
                                      num.mesh = unique(to.merge$st_mesh_nodes)))

# drop normal with NA nv and drop binom with non-NA nv
loop.params <- rbind(loop.params[(dl=='binom'  & is.na(nv)),],
                     loop.params[(dl=='normal' & !is.na(nv)),])
loop.params[dl=='binom', nv:=0.01] # default for binomial

cov.widths <- 0.8#c(0.5, 0.8, 1.0) ## set the coverage interval width
compar.dir <- "/home/merms/Dropbox/AaronJon/TMB_SpatialStats/figures/cont_sim_res/"
for(cov.width in cov.widths){
  for(ii in 1:loop.params[,.N]){
    for(with.tmb in "both"){# }c("only", "not", "both")){

      d.l <- dl  <- loop.params[ii, dl]
      n.v <- nv  <- loop.params[ii, nv]
      inla.strat <- loop.params[ii, inla.strat]
      num.mesh   <- loop.params[ii, num.mesh]

      # we only need to plot the TMB only one once, not once for every INLA!
      if(with.tmb=="only" & inla.strat != loop.params[1, inla.strat]){
        next
      }

      message(sprintf('-- plotting with.tmb = %s and %0.0f%% ints: dl=%s + n.v=%02f + inla.strat = %s + num.mesh = %i',
                      with.tmb, cov.width * 100, dl, nv, inla.strat, num.mesh))

      ## make the plot title
      if(d.l == 'normal'){
        pt <- sprintf('Average Pixel Coverage across Spatial Domain, Stratified by GP Decile, for %s Observations with Var = %.03f:\n',
                      stringr::str_to_title(d.l), n.v)
      }else{
        pt <- sprintf('Average Pixel Coverage across Spatial Domain, Stratified by GP Decile, for %s Observations:\n',
                      stringr::str_to_title(d.l))
      }
      pt <- paste(pt, sprintf('%i SPDE Vertices, Median and  %0.0f%% Quantile Bars across Monte Carlo Replicates',
                              num.mesh, cov.width*100))

      # if(d.l == 'binom') {
      #   cov.gp.sub <- subset(cov.gp, dl == d.l) ## NOTE norm.var==0.01 just to make sure nothing breaks
      # }else if(d.l == 'normal') {
      #   cov.gp.sub <- subset(cov.gp, dl == d.l & nv == n.v)
      # }

      # subset to particular liklihood and inla hyper param integration type
      if(with.tmb=="only"){
        ## only plot TMB
        cov.gp.sub <- subset(cov.gp,
                             dl == d.l &
                               nv == n.v &
                               st_mesh_nodes == num.mesh &
                               grepl("TMB", cov.gp$fit_type))
      }
      if(with.tmb=="not"){
        # only plot INLA results
        cov.gp.sub <- subset(cov.gp,
                             dl == d.l &
                               nv == n.v &
                               st_mesh_nodes == num.mesh &
                               grepl(inla.strat, cov.gp$fit_type))

      }
      if(with.tmb=="both"){
        # plot INLA and TMB results
        cov.gp.sub <- subset(cov.gp,
                             dl == d.l &
                               nv == n.v &
                               st_mesh_nodes == num.mesh &
                               (grepl(inla.strat, cov.gp$fit_type)  |
                                  grepl("TMB", cov.gp$fit_type)))
      }

      cov.gp.sub[fit_type == 'INLA_CCD_L', fit_type := 'INLA: CCD + L']
      cov.gp.sub[fit_type == 'INLA_CCD_G', fit_type := 'INLA: CCD + G']
      cov.gp.sub[fit_type == 'INLA_CCD_S', fit_type := 'INLA: CCD + S']
      cov.gp.sub[fit_type == 'INLA_EB_S', fit_type := 'INLA: EB + S']
      cov.gp.sub[fit_type == 'INLA_EB_L', fit_type := 'INLA: EB + L']
      cov.gp.sub[fit_type == 'INLA_EB_G', fit_type := 'INLA: EB + G']

      ## facet by ???
      cov.gp.sum <- summarySE(cov.gp.sub, measurevar="obs_cov",
                              groupvars=c("dl",'nv', 'n.clust', "fit_type", 'clust.var', 'n_target_cov', 'gp.dec'),
                              conf.interval = cov.width)

      ## set NA clust var to 0
      cov.gp.sum$clust.var[which(is.na(cov.gp.sum$clust.var))] <- 0

      ## groups to plot lines
      cov.gp.sum$line_group <- paste(cov.gp.sum$fit_type, cov.gp.sum$gp.dec, sep='_')

      ## rename for plotting
      setnames(cov.gp.sum, 'n.clust',   'Number of Clusters')
      setnames(cov.gp.sum, 'clust.var', 'Cluster Variance')

      pd <- position_dodge(0.05)
      gp_cov_summary <- ggplot(cov.gp.sum,
                               aes(x=n_target_cov, y=med,
                                   shape = fit_type,
                                   linetype = fit_type,
                                   colour=factor(gp.dec),
                                   group = line_group)) +
        geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
        geom_line(position=pd) +
        geom_point(position=pd, size=2) +
        scale_colour_brewer(type='div', palette = 'RdYlBu', direction=-1) + ## RdYlBu, RdYlGn, Spectral
        geom_abline(intercept = 0, slope=1) +
        facet_grid(`Number of Clusters`~`Cluster Variance`, scales='free_y', labeller=label_both) +
        ggtitle(pt) +
        ## fix the legends a bit
        labs(color = 'GP Decile', shape='Method', linetype = 'Method') + ## legend titles
        xlab(expression(Target~Coverage:~alpha)) +
        ylab('Mean Coverage of Field Estimates') +
        theme_dark()

      ## if(interactive()){ ## then we can view
      ##   print(gp_cov_summary)
      ## }

      ggsave(sprintf('%s/plot_04_%s%s%s_gp_coverage_decile_summary_%s%s%s_%i_intervals.png',
                     compar.dir,
                     dl,
                     ifelse(dl=="normal", sprintf('_%0.2fnormVar', nv), ''),
                     sprintf('_%i_meshnodes',num.mesh),
                     ifelse(with.tmb %in% c("only", "both"), "_TMB", ""),
                     ifelse(with.tmb %in% c("both"), "_and", ""),
                     ifelse(with.tmb != "only", sprintf("_INLA_%s", inla.strat), ""),
                     cov.width*100),
             plot = gp_cov_summary,
             device = 'png', units = 'in',
             width = 10, height = 12)
    }
  }
}


## ##########################################################################
## ##########################################################################
## PLOT COVERAGE BY DISTANCE TO DATA OBS
## ##########################################################################
## ##########################################################################

message('plotting coverage by dist to data')

## remove large obj for mem
rm(cov.gp); rm(cov.gp.sub); rm(cov.gp.sum)

cov.dist        <- fread(sprintf('%sall_dist_coverage_metrics.csv', compar.dir))
cov.dist[,V1:=NULL]
cov.dist <- merge(cov.dist, to.merge, by='lvid', sort=FALSE)

cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_EB_G']
cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_EB_S']
cov.dist[inla.int.strat == 'eb' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_EB_L']
cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "gaussian", fit_type := 'INLA_CCD_G']
cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "simplified.laplace", fit_type := 'INLA_CCD_S']
cov.dist[inla.int.strat == 'ccd' & model == 'inla' & inla.approx == "laplace", fit_type := 'INLA_CCD_L']
cov.dist[model == 'tmb', fit_type := 'TMB']


loop.params <- data.table(expand.grid(dl=c('binom', 'normal'),
                                      nv=c(NA, loopvars[, unique(norm.var)]),
                                      inla.strat = c('EB_G',  "EB_S",  "EB_L",
                                                     'CCD_G', "CCD_S", "CCD_L"),
                                      num.mesh = unique(long.cov.sum[["Num. spatial mesh nodes"]])))

# drop normal with NA nv and drop binom with non-NA nv
loop.params <- rbind(loop.params[(dl=='binom'  & is.na(nv)),],
                     loop.params[(dl=='normal' & !is.na(nv)),])
loop.params[dl=='binom', nv:=0.01] # default for binomial

cov.widths <- c(0.5, 0.8, 1.0) ## set the coverage interval width
for(cov.width in cov.widths){
  for(ii in 1:loop.params[,.N]){
    for(with.tmb in c("only", "not", "both")){

      d.l <- dl  <- loop.params[ii, dl]
      n.v <- nv  <- loop.params[ii, nv]
      inla.strat <- loop.params[ii, inla.strat]
      num.mesh   <- loop.params[ii, num.mesh]

      # we only need to plot the TMB only one once, not once for every INLA!
      if(with.tmb=="only" & inla.strat != loop.params[1, inla.strat]){
        next
      }

      message(sprintf('-- plotting with.tmb = %s and %0.0f%% ints: dl=%s + n.v=%02f + inla.strat = %s + num.mesh = %i',
                      with.tmb, cov.width * 100, dl, nv, inla.strat, num.mesh))

      ## make the plot title
      if(dl == 'normal'){
        pt <- sprintf('Average Pixel Coverage by Dist to Data given %s Observations with Var = %.03f',
                      stringr::str_to_title(dl), nv)
      }else{
        pt <- sprintf('Average Pixel Coverage by Dist to Data given %s Observations',
                      stringr::str_to_title(dl))
      }
      pt <- paste(pt, sprintf('\n Num of Mesh Nodes = %i || %0.0f%% Intervals of Estimates across Simulations',
                              num.mesh, cov.width*100))


      # if(dl == 'binom') {
      #   cov.dist.sub <- subset(cov.dist, dl == d.l) ## NOTE norm.var==0.01 just to make sure nothing breaks
      # }else if(dl == 'normal') {
      #   cov.dist.sub <- subset(cov.dist, dl == d.l & nv == n.v)
      # }

      # subset to particular liklihood and inla hyper param integration type
      if(with.tmb=="only"){
        ## only plot TMB
        cov.dist.sub <- subset(cov.dist,
                               dl == d.l &
                                 nv == n.v &
                                 st_mesh_nodes == num.mesh &
                                 grepl("TMB", cov.dist$fit_type))
      }
      if(with.tmb=="not"){
        # only plot INLA results
        cov.dist.sub <- subset(cov.dist,
                               dl == d.l &
                                 nv == n.v &
                                 st_mesh_nodes == num.mesh &
                                 grepl(inla.strat, cov.dist$fit_type))

      }
      if(with.tmb=="both"){
        # plot INLA and TMB results
        cov.dist.sub <- subset(cov.dist,
                               dl == d.l &
                                 nv == n.v &
                                 st_mesh_nodes == num.mesh &
                                 (grepl(inla.strat, cov.dist$fit_type)  |
                                    grepl("TMB", cov.dist$fit_type)))
      }

      ## facet by ???
      cov.dist.sum <- summarySE(cov.dist.sub, measurevar="obs_cov",
                                groupvars=c("dl",'nv', 'n.clust', "fit_type", 'clust.var', 'n_target_cov', 'obs.dist'),
                                conf.interval = cov.width)

      ## set NA clust var to 0
      cov.dist.sum$clust.var[which(is.na(cov.dist.sum$clust.var))] <- 0

      ## drop distances without enough observations
      cov.dist.sum <- subset(cov.dist.sum, N >= 25)

      ## groups to plot lines
      cov.dist.sum$line_group <- paste(cov.dist.sum$fit_type, cov.dist.sum$obs.dist, sep='_')

      ## rename for plotting
      setnames(cov.dist.sum, 'n.clust',   'Number Obs')
      setnames(cov.dist.sum, 'clust.var', 'Clust Var')

      # ## set factor order of clust.var
      # clust.var <- as.character(cov.dist.sum$clust.var)
      # clust.var[is.na(clust.var)] <- 'None'
      # cov.dist.sum$clust.var <- factor(clust.var, levels = c('None', '0.01', '0.04', '0.16'))

      pd <- position_dodge(0.05)
      gp_cov_summary <- ggplot(cov.dist.sum,
                               aes(x=n_target_cov, y=med,
                                   shape = fit_type,
                                   linetype = fit_type,
                                   colour = obs.dist,
                                   group = line_group)) +
        geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
        geom_line(position=pd) +
        geom_point(position=pd, size=2) +
        scale_colour_distiller(type='seq', palette = 10) +
        geom_abline(intercept = 0, slope=1) +
        facet_grid(`Clust Var`~`Number Obs`, scales='free_y', labeller=label_both) +
        ggtitle(pt) +
        ## fix the legends a bit
        labs(color = 'Pixels to data', shape='Method', linetype = 'Method') + ## legend titles
        xlab(expression(Target~Coverage:~alpha)) +
        ylab('Mean Coverage of Field Estimates') +
        theme_dark()

      if(interactive()){ ## then we can view
        print(gp_cov_summary)
      }

      ggsave(sprintf('%s/plot_05_%s%s%s_gp_coverage_distance_summary_%s%s%s_%i_intervals.png',
                     compar.dir,
                     dl,
                     ifelse(dl=="normal", sprintf('_%0.2fnormVar', nv), ''),
                     sprintf('_%i_meshnodes',num.mesh),
                     ifelse(with.tmb %in% c("only", "both"), "_TMB", ""),
                     ifelse(with.tmb %in% c("both"), "_and", ""),
                     ifelse(with.tmb != "only", sprintf("_INLA_%s", inla.strat), ""),
                     cov.width*100),
             plot = gp_cov_summary,
             device = 'png', units = 'in',
             width = 12, height = 8)
    }
  }
}

## ##########################################################################
## ##########################################################################
## ##########################################################################
## OLDER PLOT CODE BELOW...
## ##########################################################################
## ##########################################################################
## ##########################################################################
#
# ## ##########################################################################
# ################
# ## time plots v1 - more below ##
# ################
# ## ##########################################################################
#
# ## get the coverages calculated
# nom.cov.names <- long.cov[,sort(unique(target_cov))]
# nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100
#
# ## and assign the numeric coverage for the rows
# for(i in 1:length(nom.cov.names)){
#   long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
# }
#
# ## calculate noise to spatial ratio
# long.cov$noise_spatial_ratio <- long.cov$clust.var / long.cov$sp.var
# long.cov$noise_spatial_ratio[is.na(long.cov$noise_spatial_ratio)] <- 0
#
# ## process the total fitting and predict times
# long.cov[mean.l.model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
# long.cov[mean.l.model == 'tmb', pred.time :=  as.numeric(pred_time)]
# long.cov[mean.l.model == 'inla', fit.time :=  as.numeric(fit_time)]
# long.cov[mean.l.model == 'inla', pred.time :=  as.numeric(pred_time)]
# long.cov[, total.time :=  pred.time + fit.time]
#
# ## set cluster var NA to 0. NOTE: this is fit w/o cluster var params!
# long.cov[is.na(clust.var), clust.var:=0]
#
#
# loop.params <- data.table(dl=c('binom', rep('normal', length(loopvars[, unique(norm.var)]))),
#                           nv=c(NA, loopvars[, unique(norm.var)]))
#
# cov.widths <- c(0.5, 0.8, 1.0) ## set the coverage interval width
# for(ii in 1:loop.params[,.N]){
#
#   dl <- loop.params[ii, dl]
#   nv <- loop.params[ii, nv]
#
#   ## make the plot title
#   pt <- sprintf('Average Pixel Coverage with %s Observations | %0.0f%% E.T. Intervals',
#                 stringr::str_to_title(dl), cov.width*100)
#   if(dl == 'normal'){
#     pt <- paste(pt, sprintf('\n Gaussian SD = %.03f', nv))
#   }
#
#   if(dl == 'binom') {
#     long.cov.sub <- subset(long.cov, data.lik == dl) ## NOTE norm.var==0.01 just to make sure nothing breaks
#   }else if(dl == 'normal') {
#     long.cov.sub <- subset(long.cov, data.lik == dl & norm.var == nv)
#   }
#
#   ## facet by noise:spatial var ratio
#   long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
#                              groupvars=c('n_target_cov', 'fit_type', 'n.clust', 'data.lik',
#                                          'noise_spatial_ratio'))
#
#   ## groups to plot lines
#   long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')
#
#   ## set facet names using the labeller
#   ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
#   clust.sp.pairs <- unique(long.cov, by=c("clust.var","sp.var"))[,.(clust.var, sp.var, noise_spatial_ratio)]
#   ## sort by unique noise_spatial_ratio order
#   clust.sp.pairs <- clust.sp.pairs[order(noise_spatial_ratio),]
#   facet_labs <- as_labeller(eval(parse(text =
#                                          paste0('c(', paste0('`', sort(unique(long.cov$noise_spatial_ratio)), '`', '=', '\'',
#                                                              sprintf('(Cluster Var)/(Spatial Var): %0.2f/%0.2f = %0.2f',
#                                                                      clust.sp.pairs$clust.var,
#                                                                      clust.sp.pairs$sp.var,
#                                                                      sort(unique(long.cov$noise_spatial_ratio))), '\'', sep='', collapse=','),
#                                                 ')'))))
#
#   pd <- position_dodge(0.05)
#   fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                     aes(x=n_target_cov, y=obs_cov,
#                                         shape = fit_type,
#                                         linetype = fit_type,
#                                         colour=factor(n.clust),
#                                         group = line_group),
#                                     position=position_jitter(w=0.01, h=0.01)) +
#     geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
#     geom_line(position=pd) +
#     geom_point(position=pd) +
#     geom_abline(intercept = 0, slope = 1) +
#     facet_wrap(. ~ noise_spatial_ratio, labeller = facet_labs) +
#     ggtitle(pt) +
#     ## fix the legends a bit
#     labs(color = 'Num. Obs', shape='Method', linetype = 'Method') + ## legend titles
#     xlab(expression(Target~Coverage:~alpha)) +
#     ylab('Mean Coverage of Field Estimates')
#   if(interactive()){ ## then we can view
#     print(fit_coverage_CI_summary)
#   }
#
#   ggsave(sprintf('%s/%s_%0.2fnormVar_pixel_coverage_summary_NoiseSpatialRatio.png', compar.dir, dl, nv),
#          plot = fit_coverage_CI_summary,
#          device = 'png', units = 'in',
#          width = 12, height = 8)
# }
#
# ## ##########################################################################
# ## ##########
# ## TIME
# ## ##########
# ## ##########################################################################
#
# cm.all[model == 'tmb', fit.time := as.numeric(fit_time) + as.numeric(pt_tmb_sdreport_time)]
# cm.all[model == 'tmb', pred.time :=  as.numeric(pred_time)]
# cm.all[model == 'inla', fit.time :=  as.numeric(fit_time)]
# cm.all[model == 'inla', pred.time :=  as.numeric(pred_time)]
# cm.all[, total.time :=  pred.time + fit.time]
#
# long.cov <- data.table(gather(cm.all,
#                               operation,
#                               time_s,
#                               fit.time:total.time,
#                               factor_key = TRUE))
#
#
# ## facet by observations
# long.cov.sum <-  summarySE(long.cov, measurevar="time_s",
#                            groupvars=c("operation","fit_type", 'n.clust'))
# fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                   aes(x=n.clust, y=time_s, colour=fit_type, group = fit_type)) +
#   geom_errorbar(aes(ymin=time_s-ci, ymax=time_s+ci), width=.01) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(. ~ operation) + ggtitle(sprintf('Comparison of fit time (sec) in: %s', loopvars$data.lik[1]))
#
# if(interactive()){
#   ## then we can view
#   print(fit_coverage_CI_summary)
# }
#
# ggsave(sprintf('%s/%s_fit_time_nclust.png', compar.dir, loopvars$data.lik[1]),
#        plot = fit_coverage_CI_summary,
#        device = 'png', units = 'in',
#        width = 12, height = 12)
#
# ## facet by noise to spatial signal
# long.cov.sum <-  summarySE(long.cov, measurevar="obs_cov",
#                            groupvars=c("n_target_cov","fit_type", 'noise_spatial_ratio'))
#
# fit_coverage_CI_summary <- ggplot(long.cov.sum,
#                                   aes(x=n_target_cov, y=obs_cov, colour=fit_type, group = fit_type)) +
#   geom_errorbar(aes(ymin=obs_cov-ci, ymax=obs_cov+ci), width=.025) +
#   geom_line() +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(. ~ noise_spatial_ratio) +
#   ggtitle(sprintf('Comparison of coverage in: %s, faceted by clust.var/sp.var',
#                   loopvars$data.lik[1]))
#
# ggsave(sprintf('%s/%s_coverage_summary_noise_to_spatial_var.png', compar.dir, loopvars$data.lik[1]),
#        plot = fit_coverage_CI_summary,
#        device = 'png', units = 'in',
#        width = 12, height = 12)
#
#
# ## ## average my modeling tool
# ## if(i == 1){
# ##   cm.m.tmb  <- cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model]
# ##   cm.m.inla <- cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model]
# ## }else{
# ##   cm.m.tmb  <- rbind(cm.m.tmb,
# ##                      cm[mean.l.model == 'tmb,' lapply(.SD, mean), by=mean.l.model])
# ##   cm.m.inla <- rbind(cm.m.inla,
# ##                      cm[mean.l.model == 'inla', lapply(.SD, mean), by=mean.l.model])
# ## }


# get_jid_tid_exp_iter <- function(fn){
#   fn <- file.path(main.dir, 'logs/errors', fn)
#
#   ids <- strsplit(x=fn, '.', fixed=TRUE)
#   jid <- substr(ids[[1]][2], 2, nchar(ids[[1]][2]))
#   tid <- ids[[1]][3]
#
#   log <- readLines(fn, n=10)
#   exp.line <- grep('ON EXPERIMENT LV ID', log)
#   exp <- strsplit(log[exp.line], ': ')[[1]][2]
#   iter.line <- grep('ON SIM ITER', log)
#   iter <- strsplit(log[iter.line], ': ')[[1]][2]
#
#   return(data.table(jid=jid, tid=tid, exp=exp, iter=iter))
# }
#
# res.list <- mclapply(all.log.fns, get_jid_tid_exp_iter, mc.cores = 6)
# res.dt <- rbindlist((res.list))
