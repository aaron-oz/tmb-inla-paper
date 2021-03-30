## this script tests timing of INLA parallelization
# aoz 2021

# source("/home/merms/Documents/GitRepos/tmb_inla_comp/parallel_sim/parallel_inla_sim.R")

###########
## paths ##
###########

code.repo <- "~/Documents/GitRepos/tmb_inla_comp"
setwd(code.repo)
nga.data.dir <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects/'
outs.dir <- "~/Documents/Research/2020_tmb_v_inla/parallel_sim"
plot.dir <- "~/Dropbox/AaronJon/TMB_SpatialStats/figures/"
pardiso.lic.path <- "~/R/licenses/pardiso.lic"

##########
## load ##
##########

require(raster)
require(sp)
require(glue)
require(data.table)
require(spdep)
require(ggplot2)
require(RandomFields)
require(boot)
require(TMB)
require(INLA)
inla.prune()
# set path to pardiso solver + check if it's working
inla.setOption(pardiso.license=pardiso.lic.path)
inla.pardiso.check()
# if it isn't working, run inla.pardiso() and follow instructions

# load functions
source('./realistic_sim_utils.R')
# load gaul_list, pop_raster, simple_polygon, simple_raster, subset_shape
load(file.path(nga.data.dir, 'poly_shape.rdata'))
rm(subset_shape)

############################
## set varying sim params ##
############################

num.obs     <- c(250, 2500, 10000) # c(250, 500, 1000, 2500, 5000, 10000)
mesh.res    <- c('low', 'med', 'high') # leave out 'ultra'
inla.approx <- 'gaussian'
inla.int.strat    <- 'eb'
# check number of threads (may be more than cores)
max.threads <- parallel::detectCores()
# subset desired list to what's available
cores       <- intersect(c(1, 2, 4), 1:max.threads)

# define mesh resolutions
# R2 mesh args: largest allowed triangle edge length inner, and outer
mesh.res.dict <- list('low'   = c(0.3, 5),
                      'med'   = c(0.2, 5),
                      'high'  = c(0.15, 5),
                      'ultra' = c(0.1, 5))

# define grid of experiments
if(file.exists(file.path(outs.dir, "experiments.csv"))){
  experiments <- fread(file.path(outs.dir, "experiments.csv"))
}else{
  experiments <- data.table(expand.grid(num.obs, mesh.res, cores, inla.approx, inla.int.strat, stringsAsFactors=F))
  colnames(experiments) <- c('num.obs', 'mesh.res', 'cores', 'inla.approx', 'inla.int.strat')
  experiments[, exp := 1:.N]
  fwrite(experiments, file = file.path(outs.dir, "experiments.csv"))
}

######################
## set fixed params ##
######################

n.rep        <- 5
reg          <- 'nga'
yr.list      <- 2000
alpha        <- -1      # field int
sp.range     <- 1       # matern range
sp.var       <- .25 ^ 2 # matern var
sp.alpha     <- 2       # matern smoothness
clust.var    <- .1 ^ 2  # iid observation noise
m.clust      <- 35      # mean obs per spatial loc
sample.strat <- list(obs.loc.strat='rand',
                     urban.pop.pct=5,
                     urban.strat.pct=40)
data.lik     <- 'binom'
n.draws       <- 500

## inla approx settings
inla.int.strat <- 'ccd'
inla.approx    <- 'simplified.laplace'

## priors
# mean and sd for normal prior on fixed effects (alpha and betas)
alphaj.pri <- c(0, 3) ## N(mean, sd)
# pc.prior on clust RE precision
# (u, a) s.t. P(1/sqrt(prec) > u) = a, i.e. P(SD > u) = a
clust.prec.pri <- c(.5, .05)
# prior on spde parameters: c(a, b, c, d), where
# P(sp.range < a) = b
# P(sp.sigma > c) = d
matern.pri <- c(10, .95, 1., .05) ## a, b, c, d

# make a long data.table to store the truth and the estimates
n.exp <- nrow(experiments)
n.mth <- 1 # number of methods: inla and tmb
n.inla.tim <- 3 # fit, sample, total
n.tmb.tim  <- 7 # makeAD, Metis, Normalize, Optimize, SDrep, Sample, Total
n.tim <- n.inla.tim + n.tmb.tim
total.rows <- n.exp * n.rep * n.mth * n.tim
if(file.exists(file.path(outs.dir, "master_results.csv"))){
  master.res <- fread(file.path(outs.dir, "master_results.csv"))
}else{
  master.res <- data.table(exp = rep(1:n.exp, each = total.rows / n.exp),
                           iter = rep(rep(1:n.rep, each = n.mth * n.tim), n.exp),
                           seed = rep(-999.9, total.rows),
                           stage = rep(c("Fit", "Sample", "Total",
                                         'makeAD', 'Metis', 'Normalize', 'Optimize', 'SDrep', 'Sample', "Total"),
                                       total.rows / n.tim),
                           method = rep(c(rep("inla", n.inla.tim),
                                          rep("tmb", n.tmb.tim)), total.rows / (n.mth * n.tim)),
                           time = rep(-999.9, total.rows))
  fwrite(master.res, file = file.path(outs.dir, "master_results.csv"))
}

########################################
## prep meshes to be used across exps ##
########################################

pix.pts <- rasterToPoints(simple_raster, spatial = TRUE)
## reproject sp obj to default used in covariates
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
pix.pts <- spTransform(pix.pts, CRS(geo.prj))
pix.pts@data <- data.frame(pix.pts@data, long=coordinates(pix.pts)[,1],
                           lat=coordinates(pix.pts)[,2])
pix.pts.numeric <- as.data.frame(pix.pts@data)
## get space-locs grid to predict onto
pcoords <- xyFromCell(simple_raster, which(!is.na(values(simple_raster))))

# make and store all possible meshes and A.pred matrices
if(!file.exists(file.path(outs.dir, "meshlist.RDS"))){
  mesh.list <- list()
  for(m in mesh.res){
    ## make the triangulation
    mesh.s <- inla.mesh.2d(loc.domain = as.matrix(pix.pts.numeric[,2:3]),
                           max.e = mesh.res.dict[[m]])
    spde <- inla.spde2.pcmatern(mesh=mesh.s, alpha =2,
                                prior.range = matern.pri[1:2],
                                prior.sigma = matern.pri[3:4])
    mesh.list[[m]] <- list(mesh = mesh.s,
                           A.pred = inla.spde.make.A(mesh = mesh.s,
                                                     loc = pcoords,
                                                     group = 1),
                           spde = spde,
                           space = inla.spde.make.index("space",
                                                        n.spde = spde$n.spde,
                                                        n.group = 1))
    # print out number of REs per mesh
    print(glue('the {m} resolution mesh has {mesh.list[[m]]$mesh$n} vertices'))
  }
  saveRDS(mesh.list, file = file.path(outs.dir, "meshlist.RDS"))
}
rm(pix.pts.numeric);rm(pix.pts);rm(pcoords);rm(mesh.list)

#######################################
## make some other useful transforms ##
#######################################
## convert sp.range and sp.var into sp.kappa for rspde function
sp.kappa   <- sqrt(8) / sp.range

#####################################
# now we can start the simulations! #
#####################################

# in case we started but didn't finish
uncompleted.e.i <- unique(master.res[time == -999.9, .(exp, iter)])

for(rr in 1:nrow(uncompleted.e.i)){
  e <- uncompleted.e.i[rr, exp]
  i <- uncompleted.e.i[rr, iter]

  # reload funcs (mem issues)
  source(file.path(code.repo, 'realistic_sim_utils.R'))

  # grab the mesh.list
  mesh.list <- readRDS(file = file.path(outs.dir, "meshlist.RDS"))


  # grab the experiment settings
  n.clust <- experiments[e, num.obs]
  mesh.s  <- mesh.list[[ experiments[e, mesh.res] ]]$mesh
  A.pred  <- mesh.list[[ experiments[e, mesh.res] ]]$A.pred
  space <- mesh.list[[ experiments[e, mesh.res] ]]$space # space index for INLA
  spde  <- mesh.list[[ experiments[e, mesh.res] ]]$spde
  cores   <- experiments[e, cores]

  rm(mesh.list)

  message("\n\n\n\n\n\n\n\n\n\n")
  for(mm in 1:10){
    message(glue("\n~~~~> exp {e} ({n.exp}) + iter {i} ({n.rep}): {(e-1)*n.rep+i} of {n.exp*n.rep} <~~~~\n"))
  }
  message("\n\n\n\n\n\n\n\n\n\n")

  # (re)load the master res
  master.res <- fread(file.path(outs.dir, "master_results.csv"))
  # pull the results matrix for this run
  local.res <- master.res[exp == e & iter == i, ]
  rm(master.res)

  #######################################################################################
  ## 0, we generate a seed by hashing a string unique to this experiment and iteration ##
  #######################################################################################

  e.i <- glue("e{e}i{i}")

  hashed.seed <- strtoi(e.i, 36) %% (2 ^ 31 - 1) # modulo max seed integer

  local.res[exp == e & iter == i, seed := hashed.seed]
  set.seed(hashed.seed)

  ###################################
  ## 1, simulate and prep data obj ##
  ###################################

  sim.obj <- sim.realistic.data(reg = reg,
                                year_list = yr.list,
                                data.lik = data.lik,
                                # sd.norm = sqrt(norm.var),
                                alpha = alpha,
                                betas = NULL,
                                sp.kappa = sp.kappa,
                                sp.alpha = sp.alpha,
                                # t.rho = t.rho,
                                pixel.iid.var = NULL,
                                n.clust = n.clust,
                                m.clust = m.clust,
                                clust.re.var = clust.var,
                                # covs = covs,
                                cov_layers = list(),
                                fixed.locs = NULL,
                                fixed.gp = NULL,
                                simple_raster = simple_raster,
                                simple_polygon = simple_polygon,
                                out.dir = NULL,
                                pop_raster = pop_raster,
                                obs.loc.strat   = sample.strat[['obs.loc.strat']],
                                urban.pop.pct   = sample.strat[['urban.pop.pct']],
                                urban.strat.pct = sample.strat[['urban.strat.pct']],
                                sp.field.sim.strat = 'RF',
                                verbose = FALSE,
                                exp.iter = NULL)

  # prep data and map to mesh
  dt <- sim.obj$sim.dat; ## simulated data, lat-long, year, covs, true surface
  rm(sim.obj);rm(simple_polygon)
  dt[, id := 1:.N]
  dt[, period_id := as.numeric(as.factor(dt[, year]))]
  dt.coords <- as.matrix(dt[, .(long, lat)])
  dt.pers   <- dt[, period_id]
  nperiods  <- length(yr.list)

  A.proj <- inla.spde.make.A(mesh  = mesh.s,
                             loc   = dt.coords,
                             group = dt.pers)
  rm(dt.coords)

  #################
  ## 2, fit INLA ##
  #################

  # prep for INLA

  design_matrix <- data.frame(int = rep(1, nrow(dt)),  clust.id = 1:nrow(dt))
  stack.obs <- inla.stack(tag='est',
                          data=list(Y=dt$Y), ## response
                          A=list(A.proj,1), ## proj matrix for space, 1 for design.mat
                          effects=list(
                            space,
                            design_matrix)
                          )

  formula <- formula(Y ~ -1 + int +
                       f(clust.id, model = 'iid',
                         hyper = list(prec = list(prior='pc.prec',
                                                  param = c( clust.prec.pri[1], clust.prec.pri[2])))) +
                       f(space, model = spde))

  ## fit INLA
  pre.i.fit.t <- proc.time()[3]
  i.fit <- try(inla(formula,
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
                    control.compute=list(config = TRUE,
                                         openmp.strategy="huge"), # this forces INLA to parallelize
                    family = 'binomial',
                    num.threads = cores, #
                    Ntrials = dt$N,
                    verbose = FALSE, ## this must be false to get the logfile!
                    keep = FALSE))
  post.i.fit.t <- proc.time()[3]

  ## take draws from inla
  pre.i.draw.t <- proc.time()[3]
  i.draws <- inla.posterior.sample(n.draws, i.fit,
                                   use.improved.mean = TRUE,
                                   skew.corr = TRUE,
                                   num.threads = cores)
  post.i.draw.t <- proc.time()[3]

  # clean mem
  inla.prune()
  rm(i.fit)
  rm(i.draws)
  rm(stack.obs)
  for(g in 1:3){gc()}

  ################
  ## 3, fit TMB ##
  ################

  templ <- "model_space"
  setwd("/home/merms/Documents/GitRepos/tmb_inla_comp/parallel_sim")
  TMB::compile(paste0('./', templ,".cpp"))
  dyn.load( dynlib(templ) )

  ## make sure openmp is set to what TMB will be using
  openmp(n = cores)

  ## setup data to feed into the model
  data_full <- list(num_i = nrow(dt),  # Total number of observations
                    num_s = mesh.s$n,  # Number of vertices in SPDE mesh
                    y_i   = dt[,Y],    # Number of observed deaths in the cluster
                    n_i   = dt[,N],    # Number of exposures in the cluster
                    X_alpha  = matrix(1, nrow = nrow(dt), ncol = 1),# Covariate design matrix
                    X_betas  = matrix(1, nrow = nrow(dt), ncol = 1),# Covariate design matrix
                    M0    = spde$param.inla$M0, # SPDE sparse matrix
                    M1    = spde$param.inla$M1, # SPDE sparse matrix
                    M2    = spde$param.inla$M2, # SPDE sparse matrix
                    Aproj = A.proj,             # Projection matrix
                    options = c(0, ## if 1, run adreport
                                1, ## if 1, use priors
                                1,     ## if 1, run with intercept
                                0,     ## if 1, run with covs
                                1, ## if 1, run with cluster
                                1,           ## if 0, normal data. if 1, binom data lik
                                1), ## use normalization trick?
                    cores = cores, ## number of omp threads for TMB to use
                    flag = 1, # normalization flag. when 0, prior is returned. when 1 data is included in jnll
                    norm_prec_pri = rep(-1, 2), ## not used here - gamma on log(prec)
                    clust_prec_pri = clust.prec.pri, ## gamma on log(prec)
                    alphaj_pri = alphaj.pri, ## normal
                    ## logtau_pri = spde.theta1.pri, ## logtau: normal(mean, prec)
                    ## logkappa_pri = spde.theta1.pri ## logkappa: normal(mean, prec)
                    matern_pri = matern.pri ## c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
                    )

  ## Specify starting values for TMB parameters for GP
  tmb_params <- list(alpha = 0.0, # intercept
                     betas = 0, # not used here - cov effects
                     ## log_gauss_sigma = -2, # log(data sd) if using normal dist
                     log_gauss_prec  = -2*-2, # log(data sd)*-2 to get log(prec) if using normal dist
                     log_tau   = -1.0, # Log inverse of tau (Epsilon)
                     log_kappa = -1.0, # Matern range parameter
                     ## log_clust_sigma = -2, # log of cluster sd
                     log_clust_prec = -2*-2, # (log of cluster sd)*-2 to get log(prec)
                     clust_i = rep(0, nrow(dt)), # vector of cluster random effects
                     Epsilon_s = matrix(0, nrow=mesh.s$n, ncol=1) # GP value at obs locs
                     )

  ## make a list of things that are random effects
  rand_effs <- c('Epsilon_s')
  rand_effs <- c(rand_effs, 'clust_i')
  ## NULL out params that aren't in this run and identify extra rand effs if using them
  ADmap <- list()
  ## fixed effects coefs aren't being used
  ADmap[['betas']] <- rep(factor(NA))
  ## normal data obs variance isn't used in binom
  ADmap[['log_gauss_prec']] <- factor(NA)

  ## make the autodiff generated liklihood func & gradient
  pre.ad.t <- proc.time()[3]
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   map = ADmap,
                   hessian=TRUE,
                   DLL=templ)
  post.ad.t <- pre.sa.t <- proc.time()[3]
  runSymbolicAnalysis(obj)
  post.sa.t <- pre.norm.t <- proc.time()[3]
  ## should we use the normalization flag?
  if(data_full$options[7] == 1){
    obj <- normalize(obj, flag="flag", value = 0) ## value: Value of 'flag' that signifies to not include the data term.
  }
  post.norm.t <- pre.t.fit.t <- proc.time()[3]
  ## Run optimizer
  message('------ fitting TMB')
  opt0 <- try(do.call("nlminb",
                      list(start       =    obj$par,
                           objective   =    obj$fn,
                           gradient    =    obj$gr,
                           lower = rep(-10, length(obj$par)),
                           upper = rep( 10, length(obj$par)),
                           control     =    list(trace=1))))
  post.t.fit.t <- pre.sdr.t <- proc.time()[3]
  ## Get standard errors
  SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
  post.sdr.t <- proc.time()[3]
  ## summary(SD0, 'report')

  ## ##########
  ## PREDICT ##
  ## ##########
  message('------ making TMB predictions')

  ## now we can take draws and project to space-time raster locs
  mu <- c(SD0$par.fixed,SD0$par.random)

  ## simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }

  pre.t.draw.t <- proc.time()[3]
  L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)
  t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = n.draws)
  post.t.draw.t <- proc.time()[3]

  ##########################
  ## 4, Summarize results ##
  ##########################

  # we want:
  #  the field estimates and their coverage of the truth
  #  the param estimates and their coverage of the truth

  #~~~~~~~~~~~~~~~~~~~~~~~~
  # organize INLA results ~
  #~~~~~~~~~~~~~~~~~~~~~~~~

  inla.fit.t   <- post.i.fit.t - pre.i.fit.t
  inla.draw.t  <- post.i.draw.t - pre.i.draw.t
  inla.total.t <- inla.fit.t + inla.draw.t

  inla.times   <- c(inla.fit.t, inla.draw.t, inla.total.t)

  #~~~~~~~~~~~~~~~~~~~~~~~
  # organize TMB results ~
  #~~~~~~~~~~~~~~~~~~~~~~~

  tmb.makeAD.t <- post.ad.t - pre.ad.t
  tmb.metis.t  <- post.sa.t - pre.sa.t
  tmb.norm.t   <- post.norm.t - pre.norm.t
  tmb.fit.t    <- post.t.fit.t - pre.t.fit.t
  tmb.sdRep.t  <- post.sdr.t - pre.sdr.t
  tmb.draw.t   <- post.t.draw.t - pre.t.draw.t
  tmb.total.t <- tmb.makeAD.t + tmb.metis.t + tmb.norm.t + tmb.fit.t + tmb.sdRep.t + tmb.draw.t

  tmb.times <- c(tmb.makeAD.t, tmb.metis.t, tmb.norm.t, tmb.fit.t, tmb.sdRep.t, tmb.draw.t, tmb.total.t)

  #~~~~~~~~~~~~
  # store res ~
  #~~~~~~~~~~~~

  local.res[method == 'inla', time := inla.times]
  local.res[method == 'tmb', time := tmb.times]

  # clean up -999.9s
  local.res[local.res == -999.9] <- NA

  # slot the local into the master
  master.res <- fread(file.path(outs.dir, "master_results.csv"))
  # pull the results matrix for this run
  master.res[exp == e & iter == i, ] <- local.res

  # save to checkpoint
  fwrite(master.res, file = file.path(outs.dir, "master_results.csv"))

  ## # clean mem
  ## rm(obj)
  ## rm(opt0)
  ## rm(SD0)
  ## rm(t.draws)
  ## rm(master.res)
  ## rm(A.pred)
  ## rm(spde)
  ## rm(L)
  ## rm(mu)
  ## rm(A.proj)
  ## rm(mesh.s)
  ## rm(dt)
  ## rm(data_full)
  ## rm(sim.realistic.data)
  ## for(g in 1:3){gc()}
  ## Sys.sleep(10)

  # restart this script
  restart.r(file.to.run = "/home/merms/Documents/GitRepos/tmb_inla_comp/parallel_sim/parallel_inla_sim.R")

}   # rr: uncompleted exp/


make.plots <- TRUE
if(make.plots){

  ## process and plot
  require(ggplot2)
  # pick colors
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
  # look at them
  # pie(rep(1, 25), col = c25)

  ## load in the experiment settings, NOTE! the same dimensions were
  ## used in run1 and run2, though the INLA options varied
  experiments <- fread(file.path(outs.dir, "experiments.csv"))[, inla.approx := NULL][, inla.int.strat := NULL]

  ## load both sets of results (inla ccd+sl, and inla eb+gauss & tmb)
  master.res1 <- fread(file.path(outs.dir, 'master_results_INLA_only_CCD_SL.csv'))
  master.res1 <- na.omit(master.res1)
  master.res1[stage == 'Predict', stage := 'Prediction']
  master.res2 <- fread(file.path(outs.dir, 'master_results.csv'))
  master.res2[stage == 'Sample', stage := 'Prediction']
  ## collapse relevant TMB stages down into 'fit' and merge onto the res of master.res2
  tmb.fit <- master.res2[stage %in% c('makeAD', 'Metis', 'Normalize', 'Optimize', 'SDrep') & method == 'tmb',
                         sum(time), by = .(exp, iter)]
  setnames(tmb.fit, 'V1', 'time')
  tmb.fit <- merge(tmb.fit, unique(master.res2[, .(exp, iter, seed)]), by = c('exp', 'iter'))
  tmb.fit[, stage := 'Fit'][, method := 'tmb']
  tmb.fit <- tmb.fit[, colnames(master.res2), with = F]
  master.res2 <- rbind(master.res2, tmb.fit)
  master.res2[, inla.int.strat := 'eb'][, inla.approx := 'gaussian']
  # subset tmb to fit, predict, total
  master.res2 <- subset(master.res2, stage %in% c('Fit', 'Prediction', 'Total'))

  ## find the mean times across experiment, method, process
  mean.times1 <- master.res1[, mean(time), by = .(exp, stage, method, inla.int.strat, inla.approx)]
  mean.times2 <- master.res2[, mean(time), by = .(exp, stage, method, inla.int.strat, inla.approx)]
  mean.times <- rbind(mean.times1, mean.times2)
  ## merge on number mesh and obs info
  mean.times <- merge(mean.times, experiments, by = 'exp')

  ## spruce up names for plotting
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  mean.times[, mesh.res := firstup(mesh.res)]
  ## set factor order
  mean.times[, mesh.res := factor(mesh.res, levels = c('Low', 'Med', 'High'))]
  setnames(mean.times, c('V1', 'stage', 'mesh.res', 'method'),
           c('time', 'Operation', 'Mesh Resolution', 'Method'))
  mean.times[Method == 'inla' & inla.approx == 'simplified.laplace' & inla.int.strat == 'ccd', Method := 'INLA: CCD + S']
  mean.times[Method == 'inla' & inla.approx == 'gaussian' & inla.int.strat == 'eb', Method := 'INLA: EB + G']
  mean.times[Method == 'tmb', Method := 'TMB']

  ## make a column of the number of mesh vertices
  mean.times[, 'Number of Spatial REs' := -1]
  for(mm in mesh.res){
    mean.times[`Mesh Resolution` == firstup(mm) , `Number of Spatial REs` := mesh.list[[mm]]$mesh$n]
  }


  ## set the lines
  mean.times$line_group <- paste(mean.times$Method,   mean.times$cores)

  ## make plot
  parallel_timing_summary <- ggplot(mean.times,
                                    aes(x=log(num.obs), y=time, colour=as.factor(cores), group = line_group, linetype = Method)) +
    #  geom_errorbar(aes(ymin=l.ci, ymax=u.ci), width=.01, position=pd) +
    geom_line() + # geom_line(position=pd) +
    scale_linetype_manual(values=c("solid", "twodash", "dotted")) +
    geom_point() + # geom_point(position=pd) +
    scale_x_continuous(breaks = log(sort(mean.times$num.obs)), ## add log x-axis labels
                       labels = paste0('ln(', sort(mean.times$num.obs), ')')) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'plain', family='sans')) + ## rotate x-axis ticks
    facet_grid(`Number of Spatial REs` ~ Operation,
               labeller = labeller(.rows = label_both, .cols = label_value),
               scales='free_y') +
    ggtitle('Comparison of Parallelized Fit and Predict Times') +
    ylab("Time (seconds)") +
    xlab("Number of Clusters") +
    labs(colour="Threads", linetype = 'Method') +
    scale_colour_manual(values = c25[c(20, 2, 3, 4, 5, 7, 6)[1:3]])

  if(interactive()){
    ## then we can view
    print(parallel_timing_summary)
  }

  ggsave(sprintf('%s/parallel_times.png', plot.dir),
         plot = parallel_timing_summary,
         device = 'png', units = 'in',
         width = 9, height = 9)

}
