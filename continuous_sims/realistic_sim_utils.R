## some functions to help with realistic inla/tmb simulation comparison
## written by aoz 08mar2020

## qsub_sim: function to launch sims (and sim comparisons) on the cluster
qsub_sim_array <- function(array.ind.start = NULL, # array job start index (if launching range)
                           array.ind.stop  = NULL, # array job stop index (if launching range)
                           array.ind       = NULL, # specific array job indices to launch. overrides start/stop
                           main.dir.nm, ## head dir to store all results
                           code.path,
                           singularity = 'default',
                           singularity.opts = list(SET_OMP_THREADS = cores,
                                                   SET_MKL_THREADS = cores),
                           job.name = 'sim_array_jobs',
                           mem = '20G',
                           time = '01:00:00:00',
                           queue = 'geospatial.q',
                           priority = 0, ## default = 0, lowest priority is -1023
                           concurrent.tasks = NULL, ## max number of concurrent jobs running. if NULL, no limit
                           hold.jid = NULL, ## jobid to hold on. default no hold
                           logloc = NULL ## defaults to input/output dir in main.dir/exp.lvid
){

  ## some early checks
  if(length(unique(c(length(cov_names), length(cov_measures), length(betas)))) != 1){
    messge('cov_names, cov_measrures, and betas lengths do not match! fix and rerun')
    stop()
  }

  ## set correct project based on queue
  proj <- ifelse(queue=='geospatial.q', 'proj_geo_nodes', 'proj_geospatial')

  ## make sure we have access to J if running on non geo nodes
  node.flag <- ifelse(queue != 'geospatial.q', ' -l archive=TRUE ', '')

  ## grab the shell script we want
  shell <- '/share/code/geospatial/azimmer/lbd_core/mbg_central/share_scripts/shell_sing.sh'
  sing_image <- get_singularity(image = singularity)

  ## set the loglocation for output/error files
  if(is.null(logloc)){
    logloc <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s/logs', main.dir.nm)
  }
  error_log_dir <- paste0(logloc, '/errors/')
  output_log_dir <- paste0(logloc, '/output/')

  ## Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -N ", job.name,
                 " -e ", logloc, "/errors/",
                 " -o ", logloc, "/output/",
                 " -q ", queue,
                 " -P ", proj,
                 " -p ", priority)

  ## add on job resource requests
  qsub <- paste0(qsub,
                 ' -l m_mem_free=', mem,
                 ' -l fthread=1',
                 ' -l h_rt=', time,
                 node.flag ## either ' -l geos_node=TRUE', or ''
  )

  ## add on array job indices to run
  if(is.null(array.ind)){
    # then use start/stop indices
    qsub <- paste0(qsub,
                   ' -t ',
                   array.ind.start, ':', array.ind.stop)
  }else{
    # o.w., submit only specific jobs listed
    qsub <- paste0(qsub,
                   ' -t ',
                   paste0(array.ind, collapse=','))
  }


  ## add on concurrent task limit
  if(!is.null(concurrent.tasks)){
    qsub <- paste0(qsub,
                   ' -tc ',
                   concurrent.tasks)
  }

  ## add on stuff to launch singularity
  qsub <- qsub_sing_envs(qsub, singularity.opts,
                         sing_image)

  ## add hold on jobid flag
  if(!is.null(hold.jid)) {
    qsub <- paste(qsub,
                  "-hold_jid",
                  hold.jid)
  }

  ## append shell, and code to run
  qsub <- paste(qsub,
                shell,
                code.path)

  # ## add on all remaining arguments
  # qsub <- paste(qsub,
  #               exp.lvid, ## which row in loopvars
  #               exp.iter, ## which iteration of experiment
  #               main.dir.nm, ## which dir to load from
  #               sep = " ")

  return(qsub)
}

array_tracker <- function(init.array.jid, ## array jid
                          sims.array.jid, ## array jid
                          j.dt,           ## dt of jobs with exp, iter, jid, tid columns
                          main.dir        ## dir containing common dir with tracking csvs
){

  ## get array jobid and the range of task ids for the init job
  init.jid        <- strsplit(init.array.jid, split='[.]')[[1]][1]
  init.task.range <- strsplit(strsplit(init.array.jid, split='[.]')[[1]][2], '[:]')[[1]][1]
  init.start.id   <- strsplit(init.task.range, '[-]')[[1]][1]
  init.final.id   <- strsplit(init.task.range, '[-]')[[1]][2]

  ## get array jobid and the range of task ids for the remaining jobs
  if(!is.null(sims.array.jid)){
    sims.jid        <- strsplit(sims.array.jid, split='[.]')[[1]][1]
    sims.task.range <- strsplit(strsplit(sims.array.jid, split='[.]')[[1]][2], '[:]')[[1]][1]
    sims.start.id   <- strsplit(sims.task.range, '[-]')[[1]][1]
    sims.final.id   <- strsplit(sims.task.range, '[-]')[[1]][2]
  }else{
    j.dt <- subset(j.dt, jid==init.jid)
    sims.jid <- 1
  }

  ## query scheduler for state of the jobs
  q.ret <- system("qstat", intern = TRUE)

  ## convert to a table
  q.jobs <- do.call('rbind',
                    lapply(q.ret[-(1:2)],
                           FUN = function(x) {
                             q.ret.spl <- strsplit(x, split = ' ')[[1]];
                             data.table(jid      = q.ret.spl[3],
                                        j.status = q.ret.spl[12],
                                        tid  = q.ret.spl[56])
                           }))

  ## with array jobs, we need to parse qw separately
  if(q.jobs[,sum(j.status=='qw')>0]){
    qw.ind <- q.jobs[,which(j.status=='qw')]
    for(qq in qw.ind){
      q.ret.spl <- strsplit(q.ret[-(1:2)][qq], split = ' ')[[1]][[85]]
      qw.inds   <- strsplit(q.ret.spl, '[:]')[[1]][1]
      qw.first  <- strsplit(qw.inds, '-')[[1]][1]
      qw.last   <- strsplit(qw.inds, '-')[[1]][2]
      q.jobs.qw <- data.table(jid      = q.jobs[qq, jid],
                              j.status = 'qw',
                              tid      = seq(qw.first,
                                             qw.last)
      )
      q.jobs <- rbind(q.jobs, q.jobs.qw)
    }
  }

  ## subset relevant jobs
  ## to init.jid and sims.jid
  ## to qw with tid
  q.jobs <- subset(q.jobs, jid==init.jid | jid==sims.jid)
  q.jobs <- subset(q.jobs, tid != "")
  ## map sge.task.id to tid
  task_csv <- fread(file.path(main.dir, "task_ids.csv"))
  q.jobs[,tid := as.character(unlist(task_csv[as.numeric(q.jobs[,tid]),1]))]

  ## merge status onto jid.dt
  if(!is.null(j.dt)){
    j.dt <- merge(j.dt, q.jobs, by =c('jid', 'tid'), all.x=T, all.y=F, sort=F)
  }

  ## check the job tracking csvs
  jt.info <- check_csv_trackers(main.dir = main.dir)

  if(!is.null(j.dt)){
    ## merge on csv jt info
    j.dt <-  merge(j.dt, jt.info, by = c('exp', 'iter'), all.x=T, sort=F)

    ## determine status (not started, running, failed, completed) for each job

    ## TODO something is wrong with the logic.... jobs that don't exist are coming through as errored
    ## 2020-08-08 - I think it has to do with system r/w lag.
    ##   the job has completed but when the queen node looks at the results from the hive, the queen sees a slightly out-of-date version

    ## number of jobs in queue
    j.dt[(grepl('qw', j.status) | grepl('r', j.status)), in_q := 1]
    ## in qw means not started
    j.dt[grepl('qw', j.status)  , not_started := 1]
    ## running
    j.dt[grepl('r', j.status), running := 1]
    ## errored. script number 0 means completed
    j.dt[(is.na(j.status) & script_num != 0), errored := 1]
    ## not errored
    j.dt[is.na(errored), on_track := 1]
    ## completed
    j.dt[(is.na(j.status) & script_num==0), completed := 1]

    ## swap NA for 0 for averaging
    j.dt[is.na(j.dt)] <- 0

    ## count number of iterations in q by experiment
    exp.sum <- j.dt[, lapply(.SD, sum), by=exp, .SDcols=c('in_q')]

    ## take the mean of running status
    exp.sum <- merge(exp.sum,
                     j.dt[, lapply(.SD, FUN=function(x){round(100*mean(x, na.rm=T), 2)}),
                            by=exp,
                            .SDcols=c('on_track', 'not_started', 'running',
                                      'errored', 'completed')],
                     all.x=T, all.y=F, by='exp')

    return(list(full.tracker = j.dt,
                summ.tracker = exp.sum))
  }else{

    ## there is no job tracking, so just check for percent completion by exp
    exp.sum <- jt.info[, mean(script_num==0)*100, by=exp]
    setnames(exp.sum, 'V1', 'perc_comp')
    return(list(message=sprintf('%0.2f percent of experimets completed. %0.2f percent of iterations completed',
                                exp.sum[,mean(perc_comp==100)*100],
                                jt.info[,mean(script_num==0)]*100),
                summ.tracker=exp.sum))
  }
}

## helper function returns combos of params that never worked
## helpful to debug
list_failed_loopvars <- function(loopvars,
                                 exp.tracker){
  exp.summ <- exp.tracker$summ.tracker
  exp.full <- exp.tracker$full.tracker

  ## get all loopvar columns that vary
  unique.loopvar.cols <- loopvars[,lapply(.SD, function(x){length(unique(x))}),
                                          .SDcols = 1:ncol(loopvars)]

  ## get failed loopvars
  failed.loopvars.ind <- as.numeric(subset(exp.summ, errored > 0)$exp)
  failed.loopvars <- loopvars[failed.loopvars.ind,]
  failed.unique.loopvar.cols <- loopvars[,lapply(.SD, function(x){length(unique(x))}),
                                                 .SDcols = 1:ncol(loopvars)]

  ## out of cols that vary, here are the experiments that had something fail
  failed.combos <- unique(loopvars[failed.loopvars.ind, .SD,.SDcols=which(unique.loopvar.cols > 1)])

  ##  and here are the ones that worked
  worked.combos <- unique(loopvars[-failed.loopvars.ind, .SD,.SDcols=which(unique.loopvar.cols > 1)])

  ## find params that only showed up in failed
  always.failed <- sapply(1:ncol(failed.combos),
                          function(x){
                            setdiff(unique(failed.combos[,x, with=F]), unique(worked.combos[,x,with=F]))
                          })

  message('These loopvar columns varied and the entries shown never worked')
  print(always.failed)
}


monitor_array_jobs <- function(init.array.jid,  ## array jid
                               sims.array.jid,  ## array jid
                               j.dt,            ## dt of jobs with exp, iter, jid, tid columns
                               main.dir,        ## dir containing common dir with tracking csvs
                               pause.time = 300 ## time between checking on jobs (seconds)
){
  in_q <- 1
  while(in_q > 0){
    tracker <- array_tracker(init.array.jid=init.array.jid,
                             sims.array.jid=sims.array.jid,
                             j.dt=j.dt,
                             main.dir=main.dir)

    message(paste0('\n\n', Sys.time(), '\n\n'))
    ts <- tracker[['summ.tracker']]
    tf <- tracker[['full.tracker']]
    ## print(ts, nrow(ts))
    print(tf[, c(on_track_per=mean(on_track)*100,
                  in_q=sum(in_q),
                  running=sum(running),
                  errored=sum(errored),
                  completed=sum(completed),
                  completed_per=mean(completed)*100)],
          digits=3)
    in_q <- sum(ts$in_q)
    if(in_q > 0){ Sys.sleep(pause.time) }
  }
  return(tracker)
}

## check output dir for completion after jobs are done running
check_csv_trackers <- function(main.dir,
                               j.dt = NULL, ## if not null, merge onto csv info
                               mc.cores = num.mc.cores # load in trackers in parallel
                       ){
  ## check the csvs
  jt.dir <- paste(main.dir, 'common', 'job_tracking', sep = '/')
  jt.csvs <- list.files(path = jt.dir)
  if(length(jt.csvs) == 0) { ## no files yet,
    jt.info <- data.table(exp     = character(),
                          iter    = character(),
                          sim_loop_ct = numeric(), ## sim.loop.ct from 1_run_space_sim.R
                          script_num  = numeric())
  } else {
    jt.info <- do.call('rbind',
                       mclapply(jt.csvs,
                              FUN = function(x) {
                                csv.spl <- strsplit(x, split='[.]')[[1]][1];
                                csv.spl <- strsplit(csv.spl, split='_')[[1]];
                                csv.dat <- read.csv(paste(jt.dir, x, sep='/'));
                                data.table(exp     = csv.spl[2],
                                           iter    = csv.spl[4], ## the status is appended, use last row
                                           sim_loop_ct = csv.dat[nrow(csv.dat), 1], ## sim.loop.ct from 1_run_space_sim.R
                                           script_num  = csv.dat[nrow(csv.dat), 2]) ## script number
                              }, ## func
                       mc.cores = num.mc.cores) ## mclapply
    ) ## do.call
  } ## else: length(jt.csvs) > 0

  ## merge onto j.dt if it is passed in

}

## qsub_sim: function to launch sims (and sim comparisons) on the cluster
qsub_sim <- function(exp.lvid, ## if looping through multiple experiments - i.e. row of loopvars
                     exp.iter, ## if monte carlo iteration within an experiment
                     exp.hash, ## unique 6char string to identify all jobs belonging to the same loopvar submission
                     main.dir.nm, ## head dir to store all results
                     codepath,
                     singularity = 'default',
                     singularity_opts = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores),
                     extra_name = '',
                     mem = '20G',
                     time = '01:00:00:00',
                     queue = 'geospatial.q',
                     priority = 0, ## defaults 0, can be as low as -1023 to reduce priority
                     hold.jid = NULL, ## jobid to hold on. default no hold
                     logloc = NULL ## defaults to input/output dir in main.dir/exp.lvid/
                     ){

  ## some early checks
  if(length(unique(c(length(cov_names), length(cov_measures), length(betas)))) != 1){
    messge('cov_names, cov_measrures, and betas lengths do not match! fix and rerun')
    stop()
  }

  ## set correct project based on queue
  proj <- ifelse(queue=='geospatial.q', 'proj_geo_nodes', 'proj_geospatial')

  ## make sure we have access to J if running on non geo nodes
  node.flag <- ifelse(queue != 'geospatial.q', ' -l archive=TRUE ', '')

  ## grab the shell script we want
  shell <- '/share/code/geospatial/azimmer/lbd_core/mbg_central/share_scripts/shell_sing.sh'
  sing_image <- get_singularity(image = singularity)

  ## set the loglocation for output/error files
  if(is.null(logloc)){
    logloc <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s/%06d/logs', main.dir.nm, exp.lvid)
  }
  error_log_dir <- paste0(logloc, '/errors/')
  output_log_dir <- paste0(logloc, '/output/')

  ## Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -e ", logloc, "/errors/",
                 " -o ", logloc, "/output/",
                 " -q ", queue,
                 " -P ", proj,
                 " -p ", priority)

  ## add on job resource requests
  qsub <- paste0(qsub,
                ' -l m_mem_free=', mem,
                ' -l fthread=1',
                ' -l h_rt=', time,
                node.flag ## either ' -l geos_node=TRUE', or ''
  )

  ## add on stuff to launch singularity
  qsub <- qsub_sing_envs(qsub, singularity_opts,
                         sing_image)

  ## add hold on jobid flag
  if(!is.null(hold.jid)) {
    qsub <- paste(qsub,
                  "-hold_jid",
                  hold.jid)
  }

  ## append job name, shell, and code to run
  qsub <- paste0(qsub,
                 sprintf(" -N sim_job_%s_hash_%s_exp%06d_iter%06d",
                         extra_name, exp.hash, exp.lvid, exp.iter), ## job name
                 " ", shell, " ", codepath) ## shell and code path

  ## add on all remaining arguments
  qsub <- paste(qsub,
                exp.lvid, ## which row in loopvars
                exp.iter, ## which iteration of experiment
                main.dir.nm, ## which dir to load from
                sep = " ")

  return(qsub)
}

## used to monitor a set of jobs given amatrix of job.ids
## given a jid.dt, it tracks the progress of a fleet of jobs
## if jid.dt==NULL, then this function just queries the output folders and reports completion
track.exp.iter <- function(jid.dt, main.dir) {

  ## query the system to see all jobs
  ## first two rows are fluff
  q.ret <- system("qstat", intern = TRUE)

  ## get all jobids and their status in a nice data.table
  q.jobs <- do.call('rbind',
                    lapply(q.ret[-(1:2)],
                    FUN = function(x) { q.ret.spl <- strsplit(x, split = ' ')[[1]];
                                       data.table(jid    = q.ret.spl[3],
                                                  j.status = q.ret.spl[12])
                                     }))
  if(!is.null(jid.dt)){
    ## merge status onto jid.dt
    jid.dt <- merge(jid.dt, q.jobs, by ='jid', all.x=T, all.y=F)
  }

  ## check the csvs
  jt.dir <- paste(main.dir, 'common', 'job_tracking', sep = '/')
  jt.csvs <- list.files(path = jt.dir)
  if(length(jt.csvs) == 0) { ## no files yet,
    jt.info <- data.table(exp     = character(),
                          iter    = character(),
                          sim_loop_ct = numeric(), ## sim.loop.ct from 1_run_space_sim.R
                          script_num  = numeric())
  } else {
    jt.info <- do.call('rbind',
                       lapply(jt.csvs,
                              FUN = function(x) {
                                csv.spl <- strsplit(x, split='.')[[1]];
                                csv.spl <- strsplit(x, split='_')[[1]];
                                csv.dat <- read.csv(paste(jt.dir, x, sep='/'));
                                data.table(exp     = csv.spl[2],
                                           iter    = substr(csv.spl[4], start=1, stop=4), ## the status is appended, use last row
                                           sim_loop_ct = csv.dat[nrow(csv.dat), 1], ## sim.loop.ct from 1_run_space_sim.R
                                           script_num  = csv.dat[nrow(csv.dat), 2]) ## script number
                              } ## func
                       ) ## lapply
    ) ## do.call
  } ## else: length(jt.csvs) > 0

  if(!is.null(jid.dt)){
    ## merge on csv jt info
    jid.dt <-  merge(jid.dt, jt.info, by = c('exp', 'iter'), all.x=T)

    ## determine status (not started, running, failed, completed) for each job

    ## number of jobs in queue
    jid.dt[(grepl('qw', j.status) | grepl('r', j.status)), in_q := 1]
    ## not started means job in qw
    jid.dt[grepl('qw', j.status), not_started := 1]
    ## running
    jid.dt[grepl('r', j.status), running := 1]
    ## errored
    jid.dt[(is.na(j.status) & (script_num != 0 | is.na(script_num))), errored := 1]
    ## not errored
    jid.dt[is.na(errored), on_track := 1]
    ## completed
    jid.dt[(is.na(j.status) & script_num==0), completed := 1]

    ## swap NA for 0 for averaging
    jid.dt[is.na(jid.dt)] <- 0

    ## count number of iterations in q by experiment
    exp.sum <- jid.dt[, lapply(.SD, sum), by=exp, .SDcols=c('in_q')]

    ## take the mean of running status
    exp.sum <- merge(exp.sum,
                     jid.dt[, lapply(.SD, FUN=function(x){round(100*mean(x, na.rm=T), 2)}),
                            by=exp,
                            .SDcols=c('on_track', 'not_started', 'running',
                                      'errored', 'completed')],
                     all.x=T, all.y=F, by='exp')

    return(list(full.tracker = jid.dt,
                summ.tracker = exp.sum))
  }else{
    ## there is no job tracking, so just check for percent completion by exp
    exp.sum <- jt.info[, mean(script_num==0)*100, by=exp]
    setnames(exp.sum, 'V1', 'perc_comp')
    return(list(message=sprintf('%0.2f percent of experimets completed. %0.2f percent of iterations completed',
                                exp.sum[,mean(perc_comp==100)*100],
                                jt.info[,mean(script_num==0)]*100),
                summ.tracker=exp.sum))

  }


}

## a function to simulate from a GF using INLA
## this function is taken from the spde tutorial
rspde <- function (coords, kappa, variance = 1, alpha = 2, n = 1, mesh,
                   verbose = FALSE, seed, return.attributes = FALSE)
{

  if(is.null(seed)) seed = 0 ## If ‘seed=0L’ then GMRFLib will set the seed
                             ## intelligently/at 'random'

    t0 <- Sys.time()
    theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
    if (verbose)
        cat("theta =", theta, "\n")
    if (missing(mesh)) {
        mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords)/2)/kappa
        if (verbose)
            cat("mesh.pars =", mesh.pars, "\n")
        attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),
                                                        ],
                                               max.edge = mesh.pars[1:2],
                                               cutoff = mesh.pars[3],
                                               offset = mesh.pars[4:5]))
        if (verbose)
            cat("n.mesh =", attributes$mesh$n, "\n")
    }
    else attributes <- list(mesh = mesh)
    attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
    attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
    attributes$A <- inla.mesh.project(mesh = attributes$mesh,
                                      loc = coords)$A
    if (n == 1)
        result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                     constr = attributes$spde$f$extraconstr))
    t1 <- Sys.time()
    result <- inla.qsample(n, attributes$Q, seed = ifelse(missing(seed),
                                                          0, seed),
                           constr = attributes$spde$f$extraconstr)
    if (nrow(result) < nrow(attributes$A)) {
        result <- rbind(result, matrix(NA, nrow(attributes$A) -
                                           nrow(result), ncol(result)))
        dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
        for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*%
                                                      result[1:ncol(attributes$A),
                                                             j])
    }
    else {
        for (j in 1:ncol(result)) result[1:nrow(attributes$A),
                                         j] <- drop(attributes$A %*% result[, j])
        result <- result[1:nrow(attributes$A), ]
    }
    t2 <- Sys.time()
    attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 -
                                                              t0)
    if (return.attributes)
        attributes(result) <- c(attributes(result), attributes)
    return(drop(result))
}


sim.realistic.data <- function(reg,
                               year_list,
                               data.lik, ## either 'binom' or 'normal'
                               sd.norm = NULL, ## sd of normal observations
                               alpha = NULL, ## global intercept. if null, no intercept
                               betas = NULL, ## if null, use no covs
                               sp.kappa,
                               sp.alpha,
                               t.rho,
                               pixel.iid.var = NULL, ## pixel iid variance (spatial discontinuity)
                               n.clust,
                               m.clust,
                               clust.re.var = NULL, ## iid RE variance for each cluster observed
                               covs = NULL,
                               cov_layers = NULL, ## if supplied, use this instead of reloading covs
                               fixed.locs = NULL, ## if supplied, use these locaations to sample data
                               fixed.gp   = NULL, ## if supplied, use this GP
                               simple_raster,
                               simple_polygon,
                               out.dir,
                               pop_raster = NULL,
                               obs.loc.strat = 'rand', ## either 'rand' or 'pop.strat'. NOTE: random is proportional to population!
                               urban.pop.pct = 1, ## percent (in percent = alpha*100% space - i.e. 1 for 1%) of population that comprises urban
                               urban.strat.pct = 40, ## percent of sample locations that should come from urban pixels
                               sp.field.sim.strat = 'RF', ## one of RF or SPDE ## TODO add t-dist, extremal
                               seed = NULL,
                               verbose = FALSE,
                               exp.iter = 1){ ## exp.iter, used for file saving

  ## make some checks and set things
  if(is.null(pop_raster) & obs.loc.strat == 'pop.strat'){
    stop("You need to supply a pop raster to use obs.loc.strat=='pop.strat'")
  }
  if(obs.loc.strat == 'pop.strat' & (is.null(urban.pop.pct) | is.null(urban.strat.pct))){
      stop("You must pass in urban.pop.pct and urban.strat.pct values when you set obs.loc.strat=='pop.strat'")
  }
  if(data.lik == 'normal' & is.null(sd.norm)){
        stop("You must pass in sd.norm if setting data.lik=='normal'")
  }
  if(!is.null(betas)){
    if(!is.null(cov_layers)){
      if(length(cov_layers) != length(betas)){
        stop('The supplied cov_layers object does not match the beta arg in length. Something has gone wrong')
      }
    }else if(!is.null(covs)){
      if(nrow(covs) != length(betas)){
        stop('nrow(covs) does not match the beta arg in length. Something has gone wrong')
      }
    }else{
      stop('You must pass in either but no covs of cov_layers when you pass in betas.')
    }
  }

  ## set seed if required
  if(!is.null(seed)) set.seed(seed)

  ## create dir for simulated objects
  if(!is.null(out.dir)){
    dir.create(sprintf('%s/simulated_obj/', out.dir), recursive = T, showWarnings = F)
  }

  ## #####################################
  ## load and prepare covariate rasters ##
  ## #####################################

  ##  if(!is.null(betas)){
  if(is.null(cov_layers)){

    if(verbose) message('\n\nLOADING COVS\n')

    fixed_effects_config <- data.table(covariate = covs$name,
                                       measure = covs$meas,
                                       release = rep("2019_06_10", nrow(covs)))

    loader <- MbgStandardCovariateLoader$new(start_year = min(year_list),
                                             end_year = max(year_list),
                                             interval = 12,
                                             covariate_config = fixed_effects_config)
    cov_layers <- loader$get_covariates(simple_polygon)

    ## loop through covs, subset to years, align with simple raster, center-scale
    for(cc in 1:length(cov_layers)) {

      ## subset to years in yearlist
      if(dim(cov_layers[[cc]])[3] > 1){
        cov_layers[[cc]] <- cov_layers[[cc]][[which( min(year_list):max(year_list) %in% year_list )]]
      }

      ## align covs with simple_raster
      if(verbose) message(sprintf("On cov %i out of %i", cc, length(cov_layers)))
      cov_layers[[cc]]  <- crop(cov_layers[[cc]], extent(simple_raster))
      cov_layers[[cc]]  <- setExtent(cov_layers[[cc]], simple_raster)
      cov_layers[[cc]]  <- mask(cov_layers[[cc]], simple_raster)

      ## center-scale covs. (TODO by year or across years?)
      cov_layers[[cc]] <- (cov_layers[[cc]] - mean(values(cov_layers[[cc]]), na.rm = T)) / sd(values(cov_layers[[cc]]), na.rm = T)
    }



  }else{
    if(verbose) message('\n\nUSING PRE-SUPPLIED AND PREPPED COVS\n')
  }

  ## plot center-scaled covariates
  if(!is.null(out.dir)){
    pdf(sprintf('%s/simulated_obj/cov_plot.pdf', out.dir), width = 16, height = 16)
    for(cc in 1:length(cov_layers)){
      if(verbose) message(sprintf('Plotting covariate: %s\n', names(cov_layers)[[cc]]))
      if(dim(cov_layers[[cc]])[3] == 1){
        if(verbose) message('--plotting single synoptic map\n')
        par(mfrow = c(1, 1))
      }else{
        if(verbose) message('--plotting time series\n')
        par(mfrow = rep( ceiling(sqrt( dim(cov_layers[[cc]])[3] )), 2))
      }
      for(yy in 1:dim(cov_layers[[cc]])[3]){
        raster::plot(cov_layers[[cc]][[yy]],
                     main = paste(names(cov_layers)[cc],
                                  ifelse(dim(cov_layers[[cc]])[3] == 1,
                                         'synoptic',
                                         year_list[yy]),
                                  sep = ': '))
      }
    }
    dev.off()
  }
  ##  }else{
  ##    if(verbose) message('\n\nNO COVS\n')
  ##  }

  ## now we can simulate our true surface

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## ##########################
  ## simulate space-time gp ##
  ## ##########################

  if(verbose) message('SIMULATE SPATIAL FIELD\n')

  ## FIRST, get the pixel coords- these are useful later too

  ## convert simple raster of our region to spatialpolygonsDF
  ## these are only the NON-NA pixels!
  pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)

  ## reproject sp obj to default used in covariates
  geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  pix.pts <- spTransform(pix.pts, CRS(geo.prj))
  ## proj4string(pix.pts) ## check proj is what we wanted

  ## to simulate, we need lat-lon locs for the entire raster
  ## get coords
  pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
                             lat=coordinates(pix.pts)[,2])
  pix.pts.numeric <- as.data.frame(pix.pts@data)

  if(is.null(fixed.gp)){ ## then we simulate GP

    ## sim using SPDE mesh ## TODO export this mesh to use in fitting
    if(sp.field.sim.strat == 'SPDE'){

      ## now we can use these coords to simulate GP from rspde()
      reg.mesh <- inla.mesh.2d(boundary = inla.sp2segment(simple_polygon),
                               loc = pix.pts@data[, 2:3],
                               max.edge = c(0.1, 5),
                               offset = c(1, 5),
                               cutoff = 0.1)
      # FOR PAPER FIGURE:
      # png('~/tmb_inla_paper/spatial_matern_varying_ranges.png', width=8, height=8, units='in', res=300)
      # set.seed(413)
      # par(mfrow=c(2, 2), mai=c(.5, .5, .5, .85))
      # for(sp.range in c(.5, 1, sqrt(8), 5)){
      #   sp.kappa = sqrt(8)/sp.range

      ## get spatial fields that are GPs across space and are indep in time
      sf.iid <- rspde(coords = cbind(pix.pts.numeric[, 2], pix.pts.numeric[, 3]),
                      kappa = sp.kappa,
                      variance = sp.var,
                      alpha = sp.alpha,
                      mesh = reg.mesh,
                      n = length(year_list),
                      seed = seed)

      # FOR PAPER FIGURE:
      #   ## rasterize and plot (for testing params)
      #   sf.rast <- rasterize(x = pix.pts@data[, 2:3],
      #                        y = simple_raster_full,
      #                        field = sf.iid)
      #   if(sp.range==.5) plot(sf.rast, col=viridis(256),
      #                         main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 0.5')))
      #   if(sp.range==1) plot(sf.rast, col=viridis(256),
      #                        main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 1')))
      #   if(sp.range==sqrt(8)) plot(sf.rast, col=viridis(256),
      #                              main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = ', sqrt(8))))
      #   if(sp.range==5) plot(sf.rast, col=viridis(256),
      #                        main=expression(paste('Matern Field: ', sigma^2 == 0.25, ', ', alpha == 2, ', ', 'range = 5')))
      # }
      # dev.off()

    }else{
      reg.mesh <- NULL
    }

    ## use random fields package on simple_raster to simulate GP for spatial field
    if(sp.field.sim.strat == 'RF'){
      model <- RMmatern(nu = sp.alpha - 1, ## from INLA book
                        scale = sqrt(2 * (sp.alpha - 1)) / sp.kappa,
                        var = sp.var)

      ## sf.iid <- geostatsp::RFsimulate(model, x = simple_raster, n = length(year_list))
      sf.iid <- RFsimulate(model, x = pix.pts.numeric[, 2], y = pix.pts.numeric[, 3], n = length(year_list), spConform = FALSE)
    }

    ## simulate t dist with low DOF
    if(sp.field.sim.strat == 't'){
      stop("sp.field.sim.strat=='t' is not yet implemented")
    }

    ## simulate extremal dist
    if(sp.field.sim.strat == 'ext'){
      stop("sp.field.sim.strat=='ext' is not yet implemented")
    }

    ## ---------
    ## introduce temporal ar1 correlation at the pixel level
    if(length(year_list) > 1){ ## then, correlate gp draws
      sf.cor <- sf.iid
      for(ii in 2:ncol(sf.cor)){
        sf.cor[, ii] <- t.rho * sf.cor[, ii - 1] + sqrt(1 - t.rho ^ 2) * sf.iid[, ii]
      }


      ## convert them to rasters
      for(cc in 1:ncol(sf.iid)){
        if(cc == 1){
          sf.rast <- rasterize(x = pix.pts@data[, 2:3],
                               y = simple_raster,
                               field = sf.cor[, cc])
        }else{
          sf.rast <- addLayer(sf.rast,
                              rasterize(x = pix.pts@data[, 2:3],
                                        y = simple_raster,
                                        field = sf.cor[, cc]))
        }
      }
    }else{ ## we have a single year, no time corr needed
      sf.rast <- rasterize(x = pix.pts@data[, 2:3],
                           y = simple_raster,
                           field = sf.iid)
    }
  }else{ ## else if fixed.gp is supplied, just use that
    sf.rast <- fixed.gp
    reg.mesh <- NULL ## since this gets returned in the final list obj
  }

  ## plot gp
  if(!is.null(out.dir)){
    pdf(sprintf('%s/simulated_obj/iter%06d_st_gp_plot.pdf', out.dir, exp.iter), width = 16, height = 16)
    par(mfrow = rep( ceiling(sqrt( dim(sf.rast)[3] )), 2))
    for(yy in 1:dim(sf.rast)[3]){
      raster::plot(sf.rast[[yy]],
                   main = paste('GP',
                                year_list[yy],
                                sep = ': '))
    }
    dev.off()
  }

  ## TODO alternatively, simulate non-Gaussian field!
  ## using RandomFields, e.g.
  ## model <- RPopitz(RMexp(), alpha=2)
  ## z1 <- RFsimulate(model, seq(0,2,0.01),seq(0,2,0.01),grid=T)
  ## plot(z1, type = 'l')

  ## now we can make the iid nugget

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## simulate IID normal draws to add as nugget

  if(!is.null(pixel.iid.var)){ #} & data.lik != 'normal'){
    ## normal plus nugget means add nugget in normal draws, not to every pixel!
    ## TODO is this right to set it up and limit it this way??

    if(verbose) message('--simulate discontinuous pixel RE\n')

    ## take rnorm() draws and convert them to rasters
    for(cc in 1:length(year_list)){
      if(cc == 1){ ## initialize
        nug.rast <- rasterize(x = pix.pts@data[, 2:3],
                              y = simple_raster,
                              field = rnorm(n = nrow(pix.pts@data),
                                            mean = 0,
                                            sd = sqrt(pixel.iid.var)))
      }else{ ## add a layer
        nug.rast <- addLayer(nug.rast,
                             rasterize(x = pix.pts@data[, 2:3],
                                       y = simple_raster,
                                       field = rnorm(n = nrow(pix.pts@data),
                                                     mean = 0,
                                                     sd = sqrt(pixel.iid.var))))
      } ## else
    } ## year loop

    ## plot nugget
    pdf(sprintf('%s/simulated_obj/nugget_plot.pdf', out.dir), width = 16, height = 16)
    par(mfrow = rep( ceiling(sqrt( dim(sf.rast)[3] )), 2))
    for(yy in 1:dim(nug.rast)[3]){
      raster::plot(nug.rast[[yy]],
                   main = paste('GP',
                                year_list[yy],
                                sep = ': '))
    }
    dev.off()

  }else{
    if(verbose) message('--no spatial discontinuity RE\n')
  }


  ## now we can make the true surface

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## #####################################################
  ## make true surface by combining cov effects and gp ##
  ## #####################################################

  ## finally, we combine the gp and the covariate effects to get our surface in link (e.g. logit if binomial) space

  ## setup the global intercept
  if(is.null(alpha)) alpha <- 0

  ## start with spatial field and intercept
  true.rast <- alpha + sf.rast

  ## add on covs
  if(!is.null(betas)){
    for(cc in 1:length(cov_layers)){ ## loop though and add on coefficients*covariates to gp raster layers
      true.rast <- true.rast + betas[cc] * cov_layers[[cc]] ## should work for both stationary and time-varying
    }
  }

  ## we append the gp to the cov_layers
  ## TODO this seems like a bad idea... adjust this to keep gp out of covs
  cov_layers[['gp']] <- sf.rast

  ## and, add nugget if desired
  if(!is.null(pixel.iid.var)) true.rast <- true.rast + nug.rast

  if(!is.null(out.dir)){
    pdf(sprintf('%s/simulated_obj/iter%06d_true_surface_plot.pdf', out.dir, exp.iter), width = 16, height = 16)
    par(mfrow = rep( ceiling(sqrt( dim(true.rast)[3] )), 2))
    for(yy in 1:dim(true.rast)[3]){
      raster::plot(true.rast[[yy]],
                   main = paste('GP',
                                year_list[yy],
                                sep = ': '))
    }
    dev.off()

    saveRDS(sprintf('%s/simulated_obj/iter%06d_true_surface_raster.rds', out.dir, exp.iter), object = true.rast)
  }

  ## now the surface simulation is done and all we need to do is simulate the data

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #####################################
  ## simulate data from true surface ##
  #####################################

  if(verbose) message('SIMULATE DATA\n')


  if(is.null(fixed.locs)){ ## simulate survey locs

    ## randomly (for now) select data observation locations across time

    ## to do this, we sample, with replacement, from the lat-longs that we used to sim the GP
    if(obs.loc.strat == 'rand'){ ## select locations totally at random
      sim.rows <- sample(x = 1:nrow(pix.pts.numeric), prob = pix.pts.numeric[, 1],
                         size = n.clust * length(year_list), replace = TRUE)
    } else{ ## stratify by "urban"/"rural"
      ## given the % of population you want to be urban, find the population value cutoff
      urban_thresh <- quantile(probs = (1 - urban.pop.pct), na.omit(values(pop_raster)))

      ## make a binary urban rural raster and get the lat-longs of the pixels
      u_r_raster <- pop_raster[[1]] ## urban is 1, rural is 0
      u_r_raster[pop_raster[[1]] < urban_thresh] <- 0
      u_r_raster[pop_raster[[1]] >= urban_thresh] <- 1

      ## convert pixels to a data frame
      u_r.pts <- rasterToPoints(u_r_raster, spatial = TRUE)
      u_r.pts@data <- data.frame(u_r.pts@data, long=coordinates(u_r.pts)[,1],
                                 lat=coordinates(u_r.pts)[,2])
      u_r.pts.numeric <- as.data.frame(u_r.pts@data)

      ## sample stratified locations
      u.rows <- sample(x = which(u_r.pts.numeric[, 1] == 1), size = round(n.clust * urban.strat.pct),
                       replace = TRUE)
      r.rows <- sample(x = which(u_r.pts.numeric[, 1] == 0), size = round(n.clust * (1 - urban.strat.pct)),
                       replace = TRUE)
      sim.rows <- c(u.rows, r.rows)
    }

    ## generate a table of simulated data at the selected locations
    sim.dat <- as.data.table(pix.pts.numeric[, -1])
    sim.dat <- sim.dat[sim.rows, ]
  }else{ ## else if fixed.locs is supplied, use that
    sim.dat <- fixed.locs
  }

  ## add in years
  sim.dat[, year := rep(year_list, each = n.clust)]

  ## extract the value of the true surface at data locations
  ## this is the true value of the linear predictor! not necessarily in logit=space...
  true_lin_pred <- numeric(nrow(sim.dat))
  for(yy in unique(year_list)){
    true_lin_pred[which(sim.dat[, year] == yy)] <- raster::extract(x = true.rast[[ which(year_list %in% yy) ]],
                                                                   y = sim.dat[year == yy, .(long, lat)])
  }

  ## add in the iid cluste RE
  if(!is.null(clust.re.var)){
    if(verbose) message('-- adding in survey cluster RE')
    cluster_re <- rnorm(n=nrow(sim.dat), mean=0, sd=sqrt(clust.re.var))
  } else{
    cluster_re <- rep(0, nrow(sim.dat))
  }

  ## sim sample size of observations
  sim.dat[, N := rpois(n = nrow(sim.dat), lambda = m.clust)]

  if(data.lik == 'binom'){
    ## p_true is spatial + cluster!
    ## not the same as what the values from the truth raster
    sim.dat[, p_true := inv.logit(true_lin_pred + cluster_re)]

    ## and now we simulate binomial observations from the true surface
    sim.dat[, Y := rbinom(n = nrow(sim.dat), size = sim.dat[, N], prob = sim.dat[, p_true])]

    ## and get empirical p_obs
    sim.dat[, p_obs := Y / N]

    ## lastly, we tag on weight. TODO flesh out weighting options
    sim.dat[, weight := 1]

  }else if(data.lik == 'normal'){
    sim.dat[, p_true := true_lin_pred + cluster_re] ## identity transform

    ## and now we simulate normal observations from the true surface
    ## each individual observation has sd(obs) = sd.norm.
    ## the sd of the mean of the cluster is sd.norm/sqrt(n)
    sim.dat[, Y := apply(.SD, 1, function(x) mean(rnorm(n = x[1], mean = x[2], sd = sd.norm))),
            .SDcols = c('N', 'p_true')]

    ## and get empirical p_obs
    sim.dat[, p_obs := Y]
  }

  ## now we just finish by making some convenience objects and saving everything

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## #########################################################################################
  ## for convenience, we also extract covariate values to the same df (and the true gp val) ##
  ## ##########################################################################################

  if(verbose) message('PREPARE AND SAVE OBJECTS\n')

  cov.mat <- matrix(ncol = length(cov_layers),
                    nrow = nrow(sim.dat))
  for( cc in 1:length(cov_layers) ){
    tmp <- numeric(nrow(sim.dat))

    if(dim(cov_layers[[cc]])[3] > 1){ ## check for time-varying
      for(ll in 1:dim(cov_layers[[cc]])[3]){
        tmp[which(sim.dat[, year] == year_list[ll])] <- raster::extract(x = cov_layers[[cc]],
                                                                        y = sim.dat[year == year_list[ll], .(long, lat)],
                                                                        layer = ll)[, 1]
      }
    }else{ ## space only rasters
      tmp <- raster::extract(x = cov_layers[[cc]], y = sim.dat[, .(long, lat)])
    }

    cov.mat[, cc] <- tmp
  }
  cov.mat <- as.data.table(cov.mat)
  setnames(cov.mat, names(cov_layers))

  ## ###################################
  ## combine into a single master df ##
  ## ###################################
  sim.dat <- cbind(sim.dat, cov.mat)

  ## and drop NAs so they don't cause issues
  sim.dat <- na.omit(sim.dat)

  ## #################################
  ## save everything we might want ##
  ## #################################
  if(!is.null(out.dir)){
    saveRDS(object = sim.dat,
            file = sprintf('%s/simulated_obj/iter%06d_sim_data.rds', out.dir, exp.iter))

    saveRDS(object = cov_layers,
            file = sprintf('%s/simulated_obj/iter%06d_cov_gp_rasters.rds', out.dir, exp.iter))

    if(sp.field.sim.strat == 'SPDE' & is.null(fixed.gp)){
      saveRDS(object = reg.mesh,
              file = sprintf('%s/simulated_obj/iter%06d_region_mesh.rds', out.dir, exp.iter))
    }
  }

  #########################
  ## return a named list ##
  #########################
  return(list(sim.dat = sim.dat,
              cov.gp.rasters = cov_layers,
              true.rast = true.rast,
              mesh_s = reg.mesh ## in case we want to use the same mesh for fitting
              ))
}



#######
## Helper function for turning an xyzt table into a raster
rasterFromXYZT <- function(table,
                           z,t){
  require(data.table)
  n_periods = length(unique(table[,t]))
  table$t=table[,t]
  res=  stack(rasterFromXYZ(as.matrix(table[t==1,c('x','y',z),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      res=addLayer(res, rasterFromXYZ(as.matrix(table[t==r,c('x','y',z),with=F])))
  return(res)
}

######
## simple correlation function
my.cor <- function(x, y){
  x <- na.omit(x)
  y <- na.omit(y)
  s.x <- sum(x)
  s.x2 <- sum(x ^ 2)
  s.y <- sum(y)
  s.y2 <- sum(y ^ 2)
  s.xy <- sum(x * y)
  n <- length(x)
  return((n * s.xy - s.x * s.y) / (sqrt(n * s.x2 - s.x ^ 2) * sqrt(n * s.y2 - s.y ^ 2)))
}

#######
## continuous-rank probability score
## NOTE! since this is a normal distribution, I calcualte it on the logit scale which should be close to gaussian in our predictions
crpsNormal <- function(truth, my.est, my.var){

  sig = sqrt(my.var)
  x0 <- (truth - my.est) / sig
  res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))

  ## sign as in Held (2008)
  res <- -res

  return(res)
}

## pc.prior for gaussian precision
dPCPriPrec <- function(tau, u, a, give_log=0){
  lambda = -log(a)/u ## P( 1/sqrt(tau) > u ) = a
  logres = log(lambda/(2.0)) - (3.0/2.0)*log(tau) - lambda*1/sqrt(tau);
  if(give_log){return(logres)}else{ return(exp(logres))}
}

## summary function for long data
## utility function from the interwebs
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


###########################
###########################
## rasterize_check_coverage and build_simple_raster_pop
## added to get PR fix more quickly to deal with shapefile version


build_simple_raster_pop <- function(subset_shape,
                                    field = NULL,
                                    raking = FALSE,
                                    link_table = modeling_shapefile_version,
                                    id_raster = NULL,
                                    pop_measure = 'total',
                                    pop_release = NULL,
                                    pop_start_year = 2000,
                                    pop_end_year = 2018,
                                    shapefile_version = modeling_shapefile_version) {

  if (is.null(field)) {
    if ('GAUL_CODE' %in% names(subset_shape@data)) field <- 'GAUL_CODE'
    if ('ADM0_CODE' %in% names(subset_shape@data)) field <- 'ADM0_CODE'
  }

  if(raking) {
    field <- 'loc_id'
    # no 'loc_id' field in the link table, so we can't use it
    link_table <- NULL
  }

  ## if unspecified, get the most recent worldpop release
  if (is.null(pop_release)) {
    helper <- CovariatePathHelper$new()
    pop_rast_path  <- helper$covariate_paths(covariates = 'worldpop',
                                             measures = pop_measure,
                                             releases = pop_release)
    pop_release <- helper$newest_covariate_release(pop_rast_path)
  }

  # To ensure correct "snap" method is used, first convert to a template raster that is masked to population
  pop_rast <- brick(paste0("/snfs1/WORK/11_geospatial/01_covariates/00_MBG_STANDARD/worldpop/total/", pop_release, "/1y/worldpop_total_1y_2010_00_00.tif"))
  template_raster <- raster::crop(pop_rast, raster::extent(subset_shape), snap = "out")

  # load in the population raster given the measure and the release
  cropped_pop <- load_worldpop_covariate(template_raster,
                                         covariate = 'worldpop',
                                         pop_measure = pop_measure,
                                         pop_release = pop_release,
                                         start_year = pop_start_year,
                                         end_year = pop_end_year,
                                         interval = 12)[['worldpop']]

  ## Fix rasterize
  initial_raster <- rasterize_check_coverage(subset_shape, cropped_pop, field = field, link_table = link_table, id_raster = id_raster,
                                             shapefile_version = shapefile_version)
  if (length(subset(subset_shape, !(get(field) %in% unique(initial_raster)))) != 0) {
    rasterized_shape <-
      raster::merge(
        rasterize_check_coverage(subset(subset_shape, !(get(field) %in% unique(initial_raster))),
                                 cropped_pop,
                                 field = field,
                                 link_table = link_table,
                                 id_raster = id_raster,
                                 shapefile_version = shapefile_version),
        initial_raster)
  }
  if (length(subset(subset_shape, !(get(field) %in% unique(initial_raster)))) == 0) {
    rasterized_shape <- initial_raster
  }
  masked_pop <- raster::mask(x = cropped_pop, mask = rasterized_shape)

  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop

  return(raster_list)

}

rasterize_check_coverage <- function(shapes, template_raster, field, ..., link_table = modeling_shapefile_version, id_raster = NULL,
                                     shapefile_version = modeling_shapefile_version) {
  # backwards-compatible behavior - just call rasterize()
  if (is.null(link_table)) return(raster::rasterize(shapes, template_raster, field = field, ...))

  # Validate arguments
  is_admin_link_table <- FALSE
  if (is.data.table(link_table)) {
    is_admin_link_table <- TRUE
    # nothing to do - already a link table loaded in memory
  } else if (R.utils::isAbsolutePath(link_table)) {
    link_table <- readRDS(link_table)
  } else if (is_admin_shapefile_string(link_table)) {
    is_admin_link_table <- TRUE
    # load link table with pre-computed ownership percentages for each pixel cell
    link_table_file <- paste0(get_admin_shape_dir(link_table), "lbd_standard_link.rds")
    link_table <- readRDS(link_table_file)
  } else {
    stop("link_table argument was neither a data.table, an admin shapefile string, or an absolute path to a RDS file.")
  }

  if (! field %in% names(link_table)) {
    msg <- paste("WARNING: rasterize_check_coverage called with field", field,
                 "which is not present in link_table. Defaulting to raster::rasterize()")
    message(msg)
    return(raster::rasterize(shapes, template_raster, field = field, ...))
  }

  # aggregate link table generically for admin 0/1/2
  # Note: we need `with=FALSE` because `field` is passed as a parameter (not a hard-coded string)
  table <- link_table[,c("pixel_id", field, "area_fraction"), with = FALSE]
  if (is_admin_link_table && field != "ADM2_CODE") {
    # sum rows; area_fraction now represents the total area coverage by ADM0/1_CODE instead of ADM2_CODE
    table <- table[, .(area_fraction = sum(area_fraction)), by = c("pixel_id", field)]
  }
  # subset table so that we have 1 entry per pixel_id - the value of `field` with the maximum
  # area_fraction value for that pixel_id
  # https://stackoverflow.com/a/24558696
  pixel_owner <- table[table[, .I[which.max(area_fraction)], by = pixel_id]$V1]
  pixel_owner <- pixel_owner[order(pixel_id)]

  # Grab id_raster that shape was built against, if null defaulting to being built against world raster
  if (is.null(id_raster)) {
    # generate world raster with pixel values for `field`
    reference_pixel_owner <- suppressWarnings(empty_world_raster(shapefile_version=="2019_12_12"))
  } else {
    # Use id raster as reference, assign as NA to later fill only with pixels contained in link table
    reference_pixel_owner <- id_raster
    values(reference_pixel_owner) <- NA
  }

  # subset to only those pixels owned by a shape we're interested in
  owned_pixels <- pixel_owner[pixel_owner[[field]] %in% shapes[[field]]]
  reference_pixel_owner[owned_pixels$pixel_id] <- owned_pixels[[field]]

  result <- raster::crop(reference_pixel_owner, template_raster, snap = "near")
  if (raster::ncell(result) != raster::ncell(template_raster)) {
    message <- paste("Error in creating result raster. Should have created a raster of shape",
                     paste(dim(result), collapse=","),
                     "but instead created a raster of shape",
                     paste(dim(template_raster), collapse=","))
    stop(message)
  }
  return(result)
}


## inla_all_hyper_postprocess() ################################################
#' @title inla_all_hyper_postprocess
#' @description Postprocess multivariate normal hyperparameters to get thier
#' normal marginals
#' @param all.hyper hyper parameters from inla fit ex: inla.result$all.hyper
#' @return updated hyperparameter object
#' @rdname inla_all_hyper_postprocess

inla_all_hyper_postprocess <- function(all.hyper){
  ## postprocess all.hyper, by converting and replacing prior = 'mvnorm' into
  ## p marginals. this is for the spde-models

  len.n <- function(param, max.dim = 10000){
    len <- function(n) {
      return (n + n^2)
    }

    len.target <- length(param)
    for(n in 1:max.dim) {
      if (len(n) == len.target) {
        return(n)
      }
    }
    stop(paste("length(param) is wrong:", len.target))
  }

  get.mvnorm.marginals = function(param){
    n <- len.n(param)
    mu <- param[1:n]
    Q <- matrix(param[-(1:n)], n, n)
    Sigma <- solve((Q + t(Q))/2.)
    return(list(mean = mu, prec = 1/diag(Sigma)))
  }

  for (i in seq_along(all.hyper$random)) {
    for(j in seq_along(all.hyper$random[[i]]$hyper)) {

      if (all.hyper$random[[i]]$hyper[[j]]$prior == "mvnorm") {
        ## replace this one, and the p-following ones, with its marginals
        m <- get.mvnorm.marginals(all.hyper$random[[i]]$hyper[[j]]$param)
        for(k in 1:length(m$mean)) {
          kk <- j + k - 1
          all.hyper$random[[i]]$hyper[[kk]]$prior <- "normal"
          all.hyper$random[[i]]$hyper[[kk]]$param <- c(m$mean[k], m$prec[k])
        }
      }
    }
  }

  return (all.hyper)
}
