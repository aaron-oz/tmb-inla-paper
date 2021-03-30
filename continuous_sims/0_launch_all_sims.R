## this script can be used to launch 1_run_simulation.R in parallel on the IHME cluster
## written by aoz
## 2020MAR08
## source('/homes/azimmer/tmb_inla_comp/0_launch_all_sims.R')

## clean up
rm(list=ls())

## DO THIS!
################################################################################
## ADD A NOTE! to help identify what you were doing with this run
logging_note <-
'TRIAL 14.1: full paper submission with RandomFields simulating the truth instead of SPDE.\n with fixed cov loader'

## make a master run_date to store all these runs in a single location
main.dir.name  <- NULL ## IF NULL, run_date is made, OW uses name given
extra.job.name <- 'full_run_rf'
################################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~ SETUP ENVIRONMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## number of mc.cores for loading in csv  trackers
num.mc.cores <- 6

## is this a test run
test <- FALSE

## specify queue, project, and job requirements
q.q   <- 'geospatial.q' ## all.q ## long.q
q.m   <- '25G' ## e.g. 10G
q.t   <- '00:3:30:00' ## DD:HH:MM:SS
q.p   <- -100 ## priority: -1023 (low) - 0 (high)
q.c   <- NULL ## number concurrent running tasks allowed. NULL == no limit

#############################################
## setup the environment for singularity R ##
#############################################

## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)

## grab libraries and functions from MBG code
setwd(core_repo)
commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, 'mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(parallel)

## Now we can switch to the TMB repo
setwd(tmb_repo)
source('./realistic_sim_utils.R')
source('./qsub_utils.R')

## name overall directory with run_date if it is not named
if(is.null(main.dir.name)) main.dir.name <- make_time_stamp(TRUE)
print(main.dir.name)

## setup the main dir to store ALL experiments - i.e. each row in loopvar is an experiment
## all will be stored in main.dir, indexed by cooresponding row of loopvars
main.dir  <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s', main.dir.name)
dir.create(main.dir, recursive=T)

## write the log note
fileConn <- file(sprintf("%s/run_notes.txt", main.dir))
writeLines(logging_note, fileConn)
close(fileConn)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~ SETUP EXPERIMENTS~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################
## setup things to loop over ##
###############################

## NOTES
## this list should align with the args in 1_run_simulation.R
## NULLs can't be passed this way, so NAs are stand-ins for NULL and this gets fixed in 1_run_space_sim.R
## NAs passed in can be used to turn things off (e.g. clust.var == NULL) removes the cluster RE


## loopvars 1: region to model over
reg <- 'nga'

## loopvars: year to use (mostly for covariates)
year_list <- 2000

## loopvars 3: vector of covariates to use in model
cov_names <- c("c('access2')", ## this is just to make sure nothing breaks... this run has no covs b/c betas is set to NA
               "c('access2', 'mapincidence')")

## loopvars 4: vector of covariate measures to use in conjunction with cov_names to load covs
cov_measures <- c("c('mean')", ## this is just to make sure nothing breaks... this run has no covs b/c betas is set to NA
                  "c('mean', 'mean')")

## loopvars 5: cov effects. either NA (NO Covs), or a vector of cov effects to use wtih cov_names
betas <- c(NA, ## this is just to make sure nothing breaks... this run has no covs b/c betas is set to NA
           "c(-.25, .25)")

## loopvars 6 ## global intercept
alpha <- -1.0

## loopvars 7: spatial range as defined by INLA folks
## units are in degrees lat-long!
## Nigeria is approx 12 degrees wide and 10 degrees tall
## kappa=sqrt(8)/sp.range, so sp.range=sqrt(8) -> kappa=1 -> log(kappa)=0 (for R^2 domain)
## so, with kappa=1, 90% of the correlation drops by 2.8 degrees, or about 1/4 of the heigth/width
sp.range <- c(1, sqrt(8))

## loopvars 8: spatial nominal field variance as defiend by INLA folks
## sp.var = 1/(4*pi*kappa^2*tau^2) (for R^2 domain)
sp.var <- c(.25 ^ 2, 0.5 ^ 2)

## loopvars 9: matern smoothness = sp.alpha - 1 -> sp.alpha = matern smooth + 1 (for R^2 domain)
sp.alpha <- 2.0

## loopvars 10: cluster RE variance. NA means no effect
clust.var <-  c(NA, (c(1, 2, 4) / 10) ^ 2)

## loopvars 11: temporal auto-correlation (NOT USED IN SPACE-ONLY MODEL)
t.rho <-  0.8

## loopvars 12: R2 mesh args: largest allowed triangle edge length inner, and outer
mesh_s_params <- c("c(0.15, 5)", "c(0.2, 5)", "c(0.3, 5)")

## loopvars 13: number of clusters to simulate per year
n.clust <- c(250, 500, 1000, 2000, 4000, 8000)

## loopvars 14: mean number of individuals sim'ed per cluster using poisson(m.clust)
m.clust <- 35

## loopvars 15
## each entry must be a character string with the syntax for a 3 element R list containing:
## 1) obs.loc.strat: (either 'rand' or 'pop.strat')
## 2) urban.pop.pct:   a number between 0 and 100. the % of population that belongs to urban pixels
## 3) urban.strat.pct: a number between 0 and 100. the % of observations that should come fom urban pixels
sample.strat <- "list(obs.loc.strat='rand',
                      urban.pop.pct=5,
                      urban.strat.pct=40)"  ## random or by population for now. ## TODO cluster design

## loopvars 16: cores to use in laun
cores <- 1

## loopvars 17: number of fitted model draws to take
ndraws <- 500

## loopvars 18: mean and sd for normal prior on fixed effects (alpha and betas)
alphaj.pri <- "c(0, 3)" ## N(mean, sd)

## loopvars 19: pc.prior on clust RE precision
## (u, a) s.t. P(1/sqrt(prec) > u) = a, i.e. P(SD > u) = a
clust.prec.pri <- "c(.5, .05)"

## loopvars 20: INLA hyperparam integration strategy. can be 'eb', 'ccd', or 'grid'
inla.int.strat <- c('eb', 'ccd')

## loopvars 21: INLA marginal posterior approx strategy: can be 'gaussian', 'simplified.laplace' (default) or 'laplace'
inla.approx <- c('gaussian', 'simplified.laplace', 'laplace')

## loopvars 22: number of times to repeat an experiment (monte carlo simulations)
n.sim <- 25

## loopvars 23: data distribution: either 'binom' or 'normal'
data.lik <- c('normal', 'binom')

## loopvars 24: ONLY FOR data.lik=='normal'. variance of INDIVIDUAL normal data obs.
norm.var <- (c(1, 2, 4) / 10) ^ 2

## loopvars 25: pc.prior on normal individual level precision
## (u, a) s.t. P(1/sqrt(prec) > u) = a, i.e. P(SD > u) = a
norm.prec.pri <- "c(.5, .05)"

## loopvars 26: bias correct the mean estimates. NOTE: applies to both INLA and TMB!!
bias.correct <- c(TRUE)

## loopvars 27: perform sd correction. NOTE!! for TMB only
sd.correct <- c(TRUE)

## loopvars 28: pc.prior on spde parameters
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri <- "c(10, .95, 1., .05)" ## a, b, c, d

## loopvars 29: fix data locations
## if set to true, then we only simulate the data locations once per experiment and
## each iteration after the first will use the same locs
fix.locs <- FALSE

## loopvars 30: fix GP
## if set to true, then we only simulate the GP once per experiment and
## each iteration after the first will use the same GP
fix.gp <- FALSE

## TODO always add all vars to exand.grid()
## NOTE: I use a named list here to ensure the columns in loopvars are named
loopvars <- data.table(expand.grid(list(reg = reg, ## 1
                             year_list = year_list,
                             cov_names = cov_names,
                             cov_measures = cov_measures,
                             betas = betas, ## 5
                             alpha = alpha,
                             sp.range = sp.range,
                             sp.var = sp.var,
                             sp.alpha = sp.alpha,
                             clust.var = clust.var, ## 10
                             t.rho = t.rho,
                             mesh_s_params = mesh_s_params,
                             n.clust = n.clust,
                             m.clust = m.clust,
                             sample.strat = sample.strat, ## 15
                             cores = cores,
                             ndraws = ndraws,
                             alphaj.pri = alphaj.pri,
                             clust.prec.pri = clust.prec.pri,
                             inla.int.strat = inla.int.strat, ## 20
                             inla.approx = inla.approx,
                             n.sim = n.sim,
                             data.lik = data.lik,
                             norm.var = norm.var,
                             norm.prec.pri = norm.prec.pri, ## 25
                             bias.correct = bias.correct,
                             sd.correct = sd.correct,
                             matern.pri = matern.pri,
                             fix.locs  = fix.locs,
                             fix.gp = fix.gp ## 30
                             )))

## add on main.dir
loopvars$main.dir <- main.dir

####  drop wasteful combinations that don't need to be run ####

## drop varying norm.var with binomial
loopvars <- loopvars[!(data.lik=='binom' & norm.var > min(norm.var)),]

## drop rows with mismatches between covariate length, measures, and betas
covs.l <- sapply(as.character(loopvars$cov_names),
               function(x){length(eval(parse(text=x)))})
meas.l <- sapply(as.character(loopvars$cov_measures),
               function(x){length(eval(parse(text=x)))})
betas.l <- sapply(as.character(loopvars$betas),
                 function(x){length(eval(parse(text=x)))})
loopvars <- loopvars[covs.l == meas.l & covs.l == betas.l, ]

## drop rows with n.clust=8000, and inla approx == 'laplace'
## they tend to take too long and use too much mem
loopvars <- loopvars[!(n.clust == 8000 & inla.approx == "laplace"),]

## if testing, run only a bit
if(test){
  loopvars <- loopvars[1:100, ]
  loopvars[, n.sim:=3]
}

message(sprintf('YOU ARE ABOUT TO LAUNCH %i EXPERIMENTS', nrow(loopvars)))
message(sprintf('-- EACH WITH %i ITERATIONS', n.sim))
message(sprintf("---- THAT'S %i JOBS!", n.sim*nrow(loopvars)))


## save loopvars to this dir to reload into the parallel env
write.table(file = paste0('/ihme/scratch/users/azimmer/tmb_inla_sim/loopvars.csv'),
            x = loopvars,
            row.names = FALSE, sep=',')
## also save it inside the main.dir for the run to log what we did
write.table(file = paste0(main.dir, '/loopvars.csv'), x = loopvars,
            row.names = FALSE, sep=',')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~ RUN EXPERIMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################################################
## initialize runs: launch first (iter)ation of each (exp)eriment ##
####################################################################

## to prep covariates and other objects that are common within an experiment across iterations,
## we will launch two array jobs
## array job1 runs the first iteration for all experiments and preps objects
## the second array job will hold on completion of all jobs in array1 and will the run all other iterations
## not the most effecient, but better than launch qsubs for each experiment and iteration!

## create complete array index for all experiments and iterations within them
array.ind <- 1:(nrow(loopvars)*n.sim)

## save reference file that converts 1:Ntasks to task ids used in the parallel model
# this is needed to, eg, relaunch a specific set of failed jobs since SGE can only take start/stop task indices
fwrite(data.table(matrix(1:nrow(loopvars), ncol=1)), file = file.path(main.dir, "task_ids.csv"))

## create array_job1 for 1st iterations
qsub.string.array1 <- qsub_sim_array(array.ind.start = 1,
                                     array.ind.stop = nrow(loopvars),
                                     # array.ind = NULL,
                                     main.dir = main.dir.name,
                                     code.path = '/homes/azimmer/tmb_inla_comp/1_run_space_sim.R',
                                     singularity = 'default',
                                     singularity.opts = NULL,
                                     job.name = 'array_init',
                                     mem = q.m,
                                     time = q.t,
                                     queue = q.q,
                                     priority = q.p,
                                     concurrent.tasks = q.c,
                                     hold.jid = NULL, ## no hold
                                     logloc = NULL)   ## defaults

## launch the job and catch the message
sub.msg    <- system(qsub.string.array1, intern=TRUE)
print(sub.msg)

## extract the job id
init.array.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
init.jid       <- strsplit(init.array.jid, split = '[.]')[[1]][1]

## make a data.table to track jobs
job.dt <- data.table('exp'  = sprintf('%06d', rep(1:nrow(loopvars), n.sim)),
                     'iter' = sprintf('%06d', rep(1:n.sim, each=nrow(loopvars))),
                     'jid'  = rep("NA", nrow(loopvars)*n.sim),
                     'tid'  = rep("NA", nrow(loopvars)*n.sim))

## log the init jobs (first sim iter jobs) into the tracker dt
job.dt[iter==sprintf('%06d', 1), c('jid', 'tid') := data.table(jid=init.jid, tid=1:nrow(loopvars))]

init.tracker.list <- list(NULL) ## list to store trackers after each submit
num.launches <- 1
init.comlete <- 0

time.init.start <- Sys.time()

## monitor init jobs and relaunch failed
Sys.sleep(60) ## helps the tracker start up sometimes
init.tracker <- monitor_array_jobs(init.array.jid=init.array.jid,
                                   sims.array.jid=NULL,
                                   j.dt=job.dt,
                                   main.dir=main.dir,
                                   pause.time = 30) ## seconds

## store the init tracker list so we can print the broken ones
init.tracker.list[[num.launches]] <- init.tracker

## if any failed, show which loopvar params always failed
list_failed_loopvars(loopvars, init.tracker)

## check to see if any jobs errored. if, so, relaunch
relaunch.q.m <- q.m # to relaunch with more mem
relaunch.q.t <- q.t # to relaunch with more time ## TODO!

while( mean(init.tracker$summ.tracker[, completed]==100) < 1 & num.launches < 3 ){
  ## find errored jobs and resubmit
  message(sprintf('%.2f%% of your experiments completed successfully',
                  mean(init.tracker$summ.tracker[, completed]==100)*100))
  message(sprintf('%.2f%% of your iterations across all experiments completed successfully',
                  mean(init.tracker$full.tracker[, errored]==0)*100))
  message('These experiments had some iterations fail:')
  print(init.tracker$summ.tracker[errored > 0,])
  # message('These are the failed experiment iterations:')
  # print(init.tracker$full.tracker[errored==1, .(exp, iter, jid, tid)])

  ## get the task ids of the failed jobs
  failed.tid <- as.numeric(init.tracker$full.tracker[errored==1, tid])

  ## save reference file that converts 1:Ntasks to task ids used in the parallel model
  fwrite(matrix(failed.tid, ncol=1), file = file.path(main.dir, "task_ids.csv"))

  ## relaunch just these jobs with more 50% more mem
  relaunch.q.m <- paste0(as.numeric(gsub('G', '', relaunch.q.m))*1.5, 'G')
  relaunch.q.t <- paste("00",
                         ceiling(as.numeric(strsplit(relaunch.q.t, ":")[[1]][2])*2),
                         "00", "00", sep=":")
  qsub.string.array1.resub <- qsub_sim_array(array.ind.start = 1,
                                             array.ind.stop = length(failed.tid),
                                             main.dir = main.dir.name,
                                             code.path = '/homes/azimmer/tmb_inla_comp/1_run_space_sim.R',
                                             singularity = 'default',
                                             singularity.opts = NULL,
                                             job.name = 'array_init',
                                             mem = relaunch.q.m,
                                             time = relaunch.q.t,
                                             queue = q.q,
                                             priority = q.p,
                                             concurrent.tasks = q.c, ## concurrent tasks
                                             hold.jid = NULL, ## no hold
                                             logloc = NULL)   ## defaults

  ## launch the job and catch the message
  sub.msg <- system(qsub.string.array1.resub, intern=TRUE)
  print(sub.msg)

  ## save the job ids
  init.array.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
  init.jid       <- strsplit(init.array.jid, split = '[.]')[[1]][1]

  ## update the tracker with the init jobs
  job.dt[iter=='000001' & tid %in% failed.tid, 'jid' := init.jid]

  Sys.sleep(60) ## helps the tracker start up sometimes
  init.tracker <- monitor_array_jobs(init.array.jid=init.array.jid,
                                     sims.array.jid=NULL,
                                     j.dt=job.dt,
                                     main.dir=main.dir,
                                     pause.time = 30) ## seconds
  ## sometimes the filesystem lags, so wait and retrack jobs in case some of the
  ##   'failures' aren't really fails but just lags on viewing the tracker csvs
  Sys.sleep(180)
  init.tracker <- monitor_array_jobs(init.array.jid=init.array.jid,
                                     sims.array.jid=NULL,
                                     j.dt=job.dt,
                                     main.dir=main.dir,
                                     pause.time = 30) ## seconds

  num.launches <- num.launches + 1
  ## store the init tracker list so we can print the broken ones
  init.tracker.list[[num.launches]] <- init.tracker

  ## if any failed, show which loopvar params always failed
  list_failed_loopvars(loopvars, init.tracker)

}

if(mean(init.tracker$summ.tracker[, completed]==100) < 1){
  for(i in 1:3){message('\nSOME INIT JOBS FAILED')}
  {stop("Stopping the launch script")}
}else{ ## otherwise, they completed!
  for(i in 1:3){message('\nALL INIT JOBS SUCCESSFUL')}
}

time.init.stop <- Sys.time()
message('running time - SIMS:')
print(time.init.stop - time.init.start)

##########################################
## launch all remaining (non-init) sims ##
##########################################

if( mean(init.tracker$summ.tracker[, completed]==100) == 1 ){

  ## save reference file that converts 1:Ntasks to task ids used in the parallel model
  fwrite(data.table(matrix(array.ind[-(1:nrow(loopvars))], ncol=1)), file = file.path(main.dir, "task_ids.csv"))

  ## create array_job2 for all other sim iterations
  qsub.string.array2 <- qsub_sim_array(array.ind.start = 1,
                                       array.ind.stop  = length(array.ind[-(1:nrow(loopvars))]),
                                       array.ind = NULL,
                                       main.dir = main.dir.name,
                                       code.path = '/homes/azimmer/tmb_inla_comp/1_run_space_sim.R',
                                       singularity = 'default',
                                       singularity.opts = NULL,
                                       job.name = 'array_sims',
                                       mem = q.m,
                                       time = q.t,
                                       queue = q.q,
                                       priority = q.p,
                                       hold.jid = init.jid, ## hold this array until other array is 100% done
                                       logloc = NULL)

  ## launch the job and catch the message
  sub.msg <- system(qsub.string.array2, intern=TRUE)
  print(sub.msg)

  ## extract the job ids
  sims.array.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
  sims.jid        <- strsplit(sims.array.jid, split = '[.]')[[1]][1]

  ## update the tracker with the all the remaining jobs
  job.dt[iter != sprintf('%06d',1), c('jid', 'tid') := data.table(jid=sims.jid, tid=array.ind[-(1:nrow(loopvars))])]

  ## set the key on tid
  setkey(job.dt, 'exp', 'iter')

  sims.tracker.list <- list(NULL) ## list to store trackers after each submit
  num.launches <- 1
  sims.complete <- 0

  time.sims.start <- Sys.time()

  ## monitor sims jobs and relaunch failed
  Sys.sleep(60) ## helps the tracker start up sometimes
  ## this will print out the summary from both the init and the sims
  sims.tracker <- monitor_array_jobs(init.array.jid=init.array.jid,
                                     sims.array.jid=sims.array.jid,
                                     j.dt=job.dt,
                                     main.dir=main.dir,
                                     pause.time = 300)

  ## after all sims jobs completed, print the broken ones
  sims.tracker.list[[num.launches]] <- sims.tracker

  ## find errored jobs and resubmit
  message(sprintf('%.2f%% of your experiments completed successfully',
                  mean(sims.tracker$summ.tracker[, completed]==100)*100))
  message(sprintf('%.2f%% of your iterations across all experiments completed successfully',
                  mean(sims.tracker$full.tracker[, errored]==0)*100))
  message('These experiments had some iterations fail:')
  print(sims.tracker$summ.tracker[errored > 0,])
  # message('These are the failed experiment iterations:')
  # print(init.tracker$full.tracker[errored==1, .(exp, iter, jid, tid)])

  relaunch.q.m <- q.m
  relaunch.q.t <- q.t
  while(num.launches < 5 & sims.complete != 1){ ## relaunch broken jobs up to 2x

    if(sims.tracker$summ.tracker[errored > 0,.N] > 0){

      ## get the task ids of the failed jobs
      failed.tid <- as.numeric(sims.tracker$full.tracker[errored==1, tid])

      ## save reference file that converts 1:Ntasks to task ids used in the parallel model
      fwrite(matrix(failed.tid, ncol=1), file = file.path(main.dir, "task_ids.csv"))

      ## relaunch just these jobs with 50% more mem and more time
      relaunch.q.m <- paste0(as.numeric(gsub('G', '', relaunch.q.m))*1.5, 'G')
      relaunch.q.t <- paste("00",
                            ceiling(as.numeric(strsplit(relaunch.q.t, ":")[[1]][2])*2),
                            "00", "00", sep=":")
      qsub.string.array2.resub <- qsub_sim_array(array.ind.start = 1,
                                                 array.ind.stop = length(failed.tid),
                                                 main.dir = main.dir.name,
                                                 code.path = '/homes/azimmer/tmb_inla_comp/1_run_space_sim.R',
                                                 singularity = 'default',
                                                 singularity.opts = NULL,
                                                 job.name = 'array_sims',
                                                 mem = relaunch.q.m,
                                                 time = relaunch.q.t,
                                                 queue = q.q,
                                                 priority = q.p,
                                                 concurrent.tasks = q.c, ## concurrent tasks
                                                 hold.jid = NULL, ## no hold
                                                 logloc = NULL)   ## defaults

      ## launch the job and catch the message
      sub.msg <- system(qsub.string.array2.resub, intern=TRUE)
      print(sub.msg)

      ## save the job ids
      sims.array.jid <- strsplit(sub.msg, split = ' ')[[1]][[3]]
      sims.jid       <- strsplit(sims.array.jid, split = '[.]')[[1]][1]

      ## update the tracker with the sims jobs
      job.dt[tid %in% failed.tid, 'jid' := sims.jid]

      Sys.sleep(60) ## helps the tracker start up sometimes
      sims.tracker <- monitor_array_jobs(init.array.jid=init.array.jid,
                                         sims.array.jid=sims.array.jid,
                                         j.dt=job.dt,
                                         main.dir=main.dir,
                                         pause.time = 300) ## seconds

      num.launches <- num.launches + 1
      ## store the init tracker list so we can print the broken ones
      sims.tracker.list[[num.launches]] <- sims.tracker
    }

    sims.complete <- mean(sims.tracker$full.tracker[, errored]==0)

  } ## TODO, get the relaunch loop working

  time.sims.stop <- Sys.time()
  message('running time - SIMS:')
  print(time.sims.stop - time.sims.start)

  ## if any failed, show which loopvar params always failed
  list_failed_loopvars(loopvars, sims.tracker)

}

###############################
## PRINT AND SAVE JOB STATUS ##
###############################

# print total time
time.stop <- Sys.time()
message('running time - TOTAL:')
print(time.stop - time.init.start)

# combine init and sims 'full.trackers'
full.tracker <- rbind(init.tracker$full.tracker, sims.tracker$full.tracker)
# summarize across iterations
summ.tracker <- full.tracker[, lapply(.SD, FUN=function(x){round(100*mean(x, na.rm=T), 2)}),
                             by=exp,
                             .SDcols=c('on_track', 'not_started', 'running',
                                       'errored', 'completed')]

# print summary of job completion
message(sprintf('%.2f%% of your experiments completed successfully',
                mean(summ.tracker[, completed]==100)*100))
message(sprintf('%.2f%% of your iterations across all experiments completed successfully',
                mean(full.tracker[, errored]==0)*100))
message('These experiments had some iterations fail:')
print(summ.tracker[errored > 0,])
message('These are the failed experiment iterations:')
print(full.tracker[errored==1, .(exp, iter, jid)])

## save the environment - mostly the job ideas and qsub commands
save(list=ls(), file = sprintf('%s/completed_env.rdata', main.dir))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~ COMBINE RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## finally, we can stitch together all the experiments and simulations
source('./7_stitch_together_results_all_sims.R')

# ## print columns of loopvar that vary
# ## so we can easily see what's going on in the experiments that fail...
loopvars[init.tracker$full.tracker[errored==1, as.numeric(unique(exp))],
         !apply(loopvars, MARGIN = 2,
                 FUN=function(x){col.var <- sort(x, decreasing=F)[1] == sort(x, decreasing=T)[1]
                 if(is.na(col.var)){
                   return(TRUE)
                 }else{
                   return(col.var)}
                 }), with=F]
