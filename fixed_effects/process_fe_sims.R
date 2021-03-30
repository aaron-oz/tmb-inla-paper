# this script processes results from  run_fe_sims.R
# to make plots comparing multiple group fixed effects model fits from TMB and R-INLA
# one intercept per group - poisson
# aoz - 2020

###############
# USER INPUTS #
###############

## root dirs

code_root <- '/home/merms/Documents/GitRepos/tmb_inla_comp/fixed_effects'
outs_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/fixed_effects'
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
source(file.path(code_root, 'util_funcs.R'))

# load in the experiments for this run

# 2020_12_18 only varied num.obs with poisson
exps <- fread("experiments.csv")

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
grp.cov <- long.cov[grep('grp', par), ] # total group (this is everything when using poisson and only group means)
par.cov <- long.cov[!grep('grp', par), ]

# average over regions
grp.cov <- grp.cov[, .(obs_cov = mean(obs_cov)), by = list(exp, iter, method, n_target_cov)]

# merge on what varied in the experiment
grp.cov <- merge(grp.cov, exps[, .(exp, num.obs)], by = "exp")
par.cov <- merge(par.cov, exps[, .(exp, num.obs)], by = "exp")

#############################################
## make plots of average field coverage    ##
## - averaging over regions and iterations ##
#############################################
cov.width <- 0.8

pt <- sprintf("Average Group Coverage across 37 Fixed Effects Groups for Poisson Observations:\n Median and %0.0f%% Quantile Bars across Monte Carlo Replicates", cov.width * 100)

# get the lower and upper quantiles needed
long.cov.sum <-  summarySE(grp.cov, measurevar="obs_cov",
                           groupvars = c('n_target_cov', 'method', 'num.obs'),
                           conf.interval = cov.width)
## groups to plot lines
long.cov.sum$line_group <- paste(long.cov.sum$method, long.cov.sum$clust.var, sep='_')

## rename for plotting
setnames(long.cov.sum, 'method',  'Method')
long.cov.sum$Method <- toupper(long.cov.sum$Method)

# set factors w/ labels so we can get phi symbol in facet label
## long.cov.sum$phi <- factor(long.cov.sum$phi,
##                            levels=c("0.25","0.5","0.75", "0.9"),
##                            labels=paste0("phi: ", c(0.25, 0.5, 0.75, 0.9)))

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
  facet_grid(num.obs.f ~ ., labeller = label_parsed)

if(interactive()){ ## then we can view
  print(fit_coverage_CI_summary)
}

ggsave(file.path(outs_dir, "fixed_effects_avg_region_coverage.png"),
       plot = fit_coverage_CI_summary,
       device = 'png', units = 'in',
       width = 10, height = 12)
