# this script launches discrete simulations to compare fits from TMB and R-INLA
# aoz - 2020
# source('~/Documents/GitRepos)

###############
# USER INPUTS #
###############

test <- FALSE

## root dirs

# data_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects'
code_root <- '/home/merms/Documents/GitRepos/tmb_inla_comp/fixed_effects'
outs_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/fixed_effects'
outs_name <- ifelse(test, "test", NA) # append to outs_root for output dir. if NULL, code uses YYYY_MM_DD

## options for simulation. the factorial selection of all combos
## comprises the set of experiments, each which are run n.rep times

num.obs <- c(16, 36, 49, 100, 400) # mean number of obs per gp
grp.var <- c(1) # total var of group - ONLY USED FOR SIMULATION
# obs.var <- c(0.1, .5, 1) # var for each obs
grp.int <- -3 # intercept of groups
n.rep   <- 25 # number of reps per experiments

# set number of draws for summarization of fits
n.draws <- 500

#############
# SETUP ENV #
#############

setwd(code_root)

# set and create output dir
outs_dir <- file.path(outs_root,
                      ifelse(is.na(outs_name),
                             gsub('-', '_', as.Date(Sys.time())), outs_name))
dir.create(outs_dir, recursive = T)

# load pkgs
require(INLA)
require(TMB)
require(raster)
require(sf)
require(tmap); tmap_options(show.messages = F)
require(data.table)
require(glue)
require(scales)

# load discrete simulation functions
source(file.path(code_root, 'util_funcs.R'))

#######################################################
# make master objects for saving truths and estimates #
#######################################################

# make the set of experiments we will run
experiments <- setDT(expand.grid(num.obs, n.rep))
colnames(experiments) <- c('num.obs', 'n.rep')
experiments[, exp := 1:.N]
fwrite(experiments, file = file.path(outs_dir, "experiments.csv"))

# make a long data.table to store the truth and the estimates
n.exp <- nrow(experiments)
n.grp <- 37 # same as nigeria
n.par <- n.grp # total effects + intercept, grp var, obs var
n.mth <- 2 # number of methods: inla and tmb
total.rows <- n.exp * n.rep * n.par * n.mth

# nominal coverages we will check:
cp <- c(25, 50, 80, 90, 95)

master.res <- data.table(exp = rep(1:n.exp, each = total.rows / n.exp),
                         iter = rep(rep(1:n.rep, each = n.mth * n.par), n.exp),
                         seed = rep(-999.9, total.rows),
                         par = rep(c(glue("grp{1:n.grp}")),
                                   total.rows / n.par),
                         method = rep(rep(c("inla", "tmb"), each = n.par), total.rows / (n.mth * n.par)),
                         truth = rep(-999.9, total.rows),
                         pop = rep(-999.9, total.rows),
                         obs = rep(-999.9, total.rows),
                         est = rep(-999.9, total.rows),
                         est.sd = rep(-999.9, total.rows))
for(c in cp){
  master.res[, glue("cov{c}") := NA]
}

fwrite(master.res, file = file.path(outs_dir, "master_results.csv"))

#####################################
# now we can start the simulations! #
#####################################

for(e in 1:nrow(experiments)){

  # grab the experiment settings
  num.obs <- experiments[e, num.obs]

  # for each iter, we simulate data, fit both INLA and TMB, and then log the results
  for(i in 1:n.rep){

    message(glue("\n~~~~> exp {e} ({n.exp}) + iter {i} ({n.rep}): {(e-1)*n.rep+i} of {n.exp*n.rep} <~~~~\n"))

    # (re)load the master res
    master.res <- fread(file.path(outs_dir, "master_results.csv"))
    # pull the results matrix for this run
    local.res <- master.res[exp == e & iter == i, ]
    rm(master.res)

    ## 0, we generate a seed by hashing a string unique to this experiment and iteration

    e.i <- glue("e{e}i{i}")

    hashed.seed <- strtoi(e.i, 36) %% (2 ^ 31 - 1)# modulo max seed integer

    local.res[exp == e & iter == i, seed := hashed.seed]
    set.seed(hashed.seed)

    ######################
    ## 1, simulate data ##
    ######################

    # simulate group means
    grp.means <- rnorm(n = n.grp, mean = grp.int, sd = sqrt(grp.var))

    # true rate
    t.rate <- exp(grp.means)

    # sampled pop for each region
    s.pop <- rpois(n = n.grp, lambda = num.obs)

    # data obs
    y.obs <- rpois(n = n.grp, lambda = t.rate * s.pop)

    # save the sim obs and field
    # the obs and truth are the same for inla and tmb, so rep(..., each=2)
    local.res[grep("grp", par), truth := rep(grp.means,  2)]
    local.res[grep("grp", par), pop := rep(s.pop, each = 2)]
    local.res[grep("grp", par), obs := rep(y.obs, each = 2)]

    #################
    ## 2, fit INLA ##
    #################

    ## prep for INLA

    # make fixed effects design matrix
    dm <- as.data.frame(diag(nrow = n.grp))
    colnames(dm) <- glue("grp{1:n.grp}")
    dm$y <- y.obs
    dm$E <- s.pop

    # set prior for the intercept
  # alpha ~ N(mean, prec)
  priors  <- list(
    alpha = list(prior = 'normal', params = c(0, 1 / 25)) # mean, prec
  )

  # Set prior on $\text{logit}(\phi)$ s.t. the mixing param prior is
  # phi ~ beta(1/2,1/2)$

  # Set prior on $\text{\tau}$ s.t. the prior no scaling stdev is
  # sigma ~ N(0, prec) with sigma > 0 (truncated)

  # hyper_priors <-  list(
  #  phi  = list(prior = "logitbeta", params=c(0.5, 0.5)),
  #  prec = list(prior = "logtgaussian", params=c(0, 1/25)) # mean, prec
  # )

  # BYM2 model, no intercept, no covs
  i.formula <- as.formula(paste(c("y ~ -1", glue("grp{1:n.grp}")), collapse = " + "))

  ## fit INLA
    invisible(capture.output(i.fit <- inla(i.formula,
                                           family = "poisson", E = E,
                data = dm,
                control.fixed = list(mean.intercept = 0, prec.intercept = 1/25),
                control.predictor = list(compute=FALSE),
                control.inla = list(strategy="simplified.laplace", fast=FALSE,
                                    int.strategy = "ccd"),
                control.compute = list(config = TRUE), # allow for bias correct
                verbose=FALSE, debug=TRUE, silent=TRUE)))

    ################
    ## 3, fit TMB ##
    ################

    ## prep for TMB
    # invisible(capture.output(TMB::compile("fe.cpp")))

    # load the model
    invisible(capture.output( dyn.load( dynlib("fe") ) ))

    # the data
    t.data <- list(N = n.grp,
                   y_i = y.obs,
                   e_i = s.pop,
                   X = diag(nrow = n.grp),
                   options = list(adreport_on = 1,
                                  normal_trick = 0,
                                  verbose = 0),
                   flag = 1, ## only used if options[['normal_trick']]==1
                   coef_pri = priors[['alpha']][['params']])# ,
    # bym2_phi_pri = hyper_priors[['phi']][['params']],
    # bym2_tau_pri = hyper_priors[['prec']][['params']])

    # starting vals
    t.params <- list(betas = rep(0, n.grp)) # int + grps
    #logit_phi = 0,
    #log_tau = 0,
    #Epsilon_s = matrix(0, ncol = 1, nrow = 2 * nrow(t.data[['Q_struc']])))

    # rand effs
    t.rand <- NULL
    # t.rand <- c('Epsilon_s')

    ## and lastly, we can NULL out params that are in the c++ template but
    ## which won't be used in this run. this allows us to build up a more
    ## complicated template that can take different options. named items
    ## in the list which as set to factor(NA) will be left out
  ADmap <- list()

  ## fit in TMB

    # make AD functional
    invisible(capture.output(
      obj <- MakeADFun(data       = t.data,
                       parameters = t.params,
                       random     = t.rand,
                       map        = ADmap,
                       hessian    = TRUE,
                       DLL        = "fe") ))

    # call the optimizer
    invisible(capture.output(
      opt0 <- try(do.call("nlminb",
                          list(start     = obj$par,
                               objective = obj$fn,
                               gradient  = obj$gr,
                               lower     = rep(-20, length(obj$par)),
                               upper     = rep( 20, length(obj$par)),
                               control   = list(trace=1)))) ))

    # report the estimates
    invisible(capture.output(
      SD0 <- TMB::sdreport(obj) ))
    # getJointPrecision=FALSE, # no REs
    # getReportCovariance = FALSE, # no REs
    # bias.correct = TRUE,
    # bias.correct.control = list(sd = TRUE))

    ## summary(SD0, report = 'report')
    #R0 <- obj$report()

    ##########################
    ## 4, Summarize results ##
    ##########################

    ## we want:
    ## the group estimates and their coverage of the truth
    ## the param estimates and their coverage of the truth

    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # organize INLA estimates ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~

    # take INLA draws

    invisible(capture.output(i.draws <-
                               suppressMessages(suppressWarnings(inla.posterior.sample(n.draws, i.fit,
                                                                                       use.improved.mean = FALSE,
                                                                                       verbose = F)))))

    # estimates of group
    i.grp.idx <- tail(grep("grp", rownames(i.draws[[1]]$latent)), n.grp)
    e.i.grp   <- data.table(do.call("cbind",
                                    lapply(i.draws,
                                           function(x){x$latent[i.grp.idx]})) )
    colnames(e.i.grp) <- glue("draw{1:n.draws}")
    e.i.grp[, par := as.character(glue("grp{1:n.grp}"))][, method := "inla"]

    # all estimates from INLA
    e.i <- rbind(e.i.grp)

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # organize TMB estimates ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~

    ## take TMB draws

    mu <- c(SD0$par.fixed)

    # first, try to get the cholesky decomp of the precision
    L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)

    #if(class(L) != "try-error"){

    # take the draws
    t.draws <- do.call("rbind", lapply(1:n.grp,
                                       FUN =  function(x){
                                         rnorm(n = n.draws, mean = mu[x], sd = sqrt(diag(SD0$cov.fixed))[x])
                                       }))

    ## separate out the tmb_draws
    t.parnames <- c(names(SD0$par.fixed), names(SD0$par.random))

    t.grp.idx <- grep("betas", t.parnames)
    e.t.grp   <- data.table(t.draws[t.grp.idx, ])
    colnames(e.t.grp) <- glue("draw{1:n.draws}")
    e.t.grp[, par := as.character(glue("grp{1:n.grp}"))][, method := "tmb"]

    # all estimates from TMB
    e.t <- rbind(e.t.grp)

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # summarize INLA and TMB ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~

    est <- rbind(e.i, e.t)

    # get est :=mean (draws)
    # get est.sd := sd(draws)
    # get covXY := quantiles(draws)

    e.sum <- est[, .(est    = apply(.SD, 1, median), # median for stability over mean
                     est.sd = apply(.SD, 1, sd)),
                 by = list(par, method),
                 .SDcols = grep("draw", names(est), value = T)]

    # find the lower and upper quantiles
    e.lui <- t(apply(est[, grep("draw", names(est), value = T), with = F], 1, quantile,
                     p = c((1 - cp/100)/2,
                           cp/100 + (1 - cp/100) / 2), na.rm=T))

    ## add the binary coverage into the res table
    for(cc in 1:length(cp)){

      ## cc is index of lower
      ## u.cc is index of upper
      u.cc <- length(cp) + cc

      li <- e.lui[, cc]
      ui <- e.lui[, u.cc]

      c <- cp[cc]

      e.sum[, glue("cov{c}") := (local.res[, truth] >= li & local.res[, truth] <= ui)]

    }

    #~~~~~~~~~~~~
    # store res ~
    #~~~~~~~~~~~~

    # slot the summary into the local resale
    slot.cols <- c("est", "est.sd", glue("cov{cp}"))
    local.res[, (slot.cols) := (e.sum[, slot.cols, with = F])]

    # clean up -999.9s
    local.res[local.res == -999.9] <- NA

    # slot the local into the master
    master.res <- fread(file.path(outs_dir, "master_results.csv"))
    # pull the results matrix for this run
    master.res[exp == e & iter == i, ] <- local.res

    # save to checkpoint
    fwrite(master.res, file = file.path(outs_dir, "master_results.csv"))

  } # i
}   # e
