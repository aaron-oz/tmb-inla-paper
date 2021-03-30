# this script processes results from  run_discrete_sims.R
# to make plots comparing fits from TMB and R-INLA
# aoz - 2020

###############
# USER INPUTS #
###############

## root dirs

data_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects'
code_root <- '/home/merms/Documents/GitRepos/tmb_inla_comp/'
outs_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/discrete_sim'
outs_name <- '2020_12_19' # append to outs_root for output dir. if NULL, code uses YYYY_MM_DD

#############
# SETUP ENV #
#############

# set and create output dir
outs_dir <- file.path(outs_root,
                      ifelse(is.null(outs_name),
                             gsub('-', '_', as.Date(Sys.time())), outs_name))
setwd(outs_dir)

# load pkgs
require(INLA)
require(TMB)
require(raster)
require(sf)
require(tmap); tmap_options(show.messages = F)
require(data.table)
require(glue)
require(scales)
require(ggplot2)
require(tidyr)
require(plyr)

# load discrete simulation functions
source(file.path(code_root, 'discrete_sim/discrete_sim_funcs.R'))

# load prepped nigeria data
load(file.path(data_root, "prepped_discrete_obj.RData"))

# convert to sf obj - used for plotting
nga1 <- st_as_sf(nga1)

# load in the experiments for this run
exps <- fread("experiments.csv")
setnames(exps, "bym.phi", "phi")
# 2020_12_03 only varied num.obs and phi

# load results from sims
res <- fread("master_results.csv")

long.cov <- data.table(gather(res,
                              target_cov,  ## name of NEW key col
                              obs_cov,     ## name of NEW value col
                              cov25:cov95, ## cols to convert from wide to long
                              factor_key = TRUE
                              ))

# get the coverages calculated
nom.cov.names <- long.cov[,sort(unique(target_cov))]
nom.cov.vals <- as.numeric(gsub("[^0-9]", "",  nom.cov.names))/100

# and assign the numeric coverage for the rows
for(i in 1:length(nom.cov.names)){
  long.cov[target_cov == nom.cov.names[i], n_target_cov := nom.cov.vals[i]]
}

# split into bym total and other params so we can average over bym
total.cov <- long.cov[grep("bym", par), ] # total bym2 coverage
param.cov <- long.cov[!grep("bym", par), ]

# average over regions
total.cov <- total.cov[, .(obs_cov = mean(obs_cov)), by = list(exp, iter, method, n_target_cov)]

# merge on what varied in the experiment
total.cov <- merge(total.cov, exps[, .(num.obs, phi, exp)], by = "exp")
param.cov <- merge(param.cov, exps[, .(num.obs, phi, exp)], by = "exp")

#############################################
## make plots of average field coverage    ##
## - averaging over regions and iterations ##
#############################################
cov.width <- 0.8

pt <- sprintf("Average Field Coverage across Spatial Regions for Poisson Observations:\n Median and %0.0f%% Quantile Bars across Monte Carlo Replicates", cov.width * 100)

# get the lower and upper quantiles needed
long.cov.sum <-  summarySE(total.cov, measurevar="obs_cov",
                           groupvars = c('n_target_cov', 'method', 'num.obs', 'phi'),
                           conf.interval = cov.width)
## groups to plot lines
long.cov.sum$line_group <- paste(long.cov.sum$method, long.cov.sum$clust.var, sep='_')

## rename for plotting
setnames(long.cov.sum, 'method',  'Method')
long.cov.sum$Method <- toupper(long.cov.sum$Method)

# set factors w/ labels so we can get phi symbol in facet label
long.cov.sum$phi <- factor(long.cov.sum$phi,
                           levels=c("0.25","0.5","0.75", "0.9"),
                           labels=paste0("phi: ", c(0.25, 0.5, 0.75, 0.9)))

long.cov.sum$num.obs.f  <- factor(long.cov.sum$num.obs,
                                  levels=c(4, 6, 7, 10, 20) ^ 2,
                                  labels=paste0("Num.~Observations:", c(4, 6, 7, 10, 20) ^ 2))

pd <- position_dodge(0.05)
fit_coverage_CI_summary <- ggplot(long.cov.sum,
                                  aes(x = n_target_cov, y = obs_cov,
                                      shape = Method,
                                      color = Method,
                                      group = Method),
                                  position=position_jitter(w=0.02, h=0.02)) +
  geom_errorbar(position=pd, aes(ymin = l.ci, ymax = u.ci), width = .025) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle(pt) +
  theme(panel.spacing.x = unit(4, "mm")) +
  ## fix the legends a bit
  ## labs(color = "Method", shape='Fit Type', linetype='Fit Type') + e## legend titles
  xlab(expression(Target~Coverage:~alpha)) +
  ylab('Mean Coverage of Field Estimates') +
  facet_grid(num.obs.f ~ phi, labeller = label_parsed)

if(interactive()){ ## then we can view
  print(fit_coverage_CI_summary)
}

ggsave(file.path(outs_dir, "avg_region_coverage.png"),
       plot = fit_coverage_CI_summary,
       device = 'png', units = 'in',
       width = 10, height = 12)

##################################
## make plots of param coverage ##
## - averaging over iterations  ##
##################################
cov.width = 0.8

param.cov[, bias := est - truth]
# param.long <- melt(param.cov[, .(num.obs, phi, par, method, est, est.sd)], measure = "par")

## make the plot title
pt <- sprintf('Parameter Bias for Poisson Observations:')
pt <- paste(pt, sprintf('Median and %0.0f%% Quantile Bars across Monte Carlo Replicates', cov.width*100))

param.long.sum <-  summarySE(param.cov, measurevar="bias",
                             groupvars = c('par', 'method', 'num.obs', 'phi'),
                             conf.interval = cov.width)

## groups to plot lines
## long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$n.clust, sep='_')

## rename for plotting
setnames(param.long.sum, 'num.obs',  'Number Obs')
setnames(param.long.sum, 'method',  'Method')
param.long.sum$Method <- toupper(param.long.sum$Method)

# set factors w/ labels so we can get phi symbol in facet label
param.long.sum$phi <- factor(param.long.sum$phi,
                             levels=c("0.25","0.5","0.75", "0.9"),
                           labels=paste0("phi: ", c(0.25, 0.5, 0.75, 0.9)))

param.long.sum$par <- factor(param.long.sum$par,
                             levels=c("int","phi","tau"),
                             labels=c("Intercept", "phi", "tau"))


pd <- position_dodge(.15)
fe_bias_summary <- ggplot(param.long.sum,
                          aes(x=log(`Number Obs`), y=med,
                              shape = Method,
                              colour = Method)) +
  geom_errorbar(position=pd, aes(ymin=l.ci, ymax=u.ci), width=.025) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2) +
  geom_abline(intercept = 0, slope=0) +
  ## facet_wrap(. ~ variable, scales='free_y') +
  ggtitle(pt) +
  ## fix the legends a bit
  ## labs(color = 'Clust. Var', shape='Fit Type', linetype = 'Fit Type') + ## legend titles
  xlab('Number of Observations per Region') +
  scale_x_continuous(breaks = log(sort(param.long.sum$`Number Obs`)), ## add log x-axis labels
                     labels = paste0('ln(', sort(param.long.sum$`Number Obs`), ')')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', family='sans')) + ## rotate x-axis ticks
  ylab('Bias: (Estimate - True)') +
  facet_grid(par ~ phi, labeller =  label_parsed, scales = "free")

if(interactive()){ ## then we can view
  print(fe_bias_summary)
}

ggsave(file.path(outs_dir, "avg_param_bias.png"),
       plot = fe_bias_summary,
       device = 'png', units = 'in',
       width = 10, height = 12)
