# this script launches discrete simulations to compare fits from TMB and R-INLA
# aoz - 2020

###############
# USER INPUTS #
###############

test <- FALSE

## root dirs

# data_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects'
code_root <- '/home/merms/Documents/GitRepos/tmb_inla_comp/re_iid'
outs_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/re_iid'
outs_name <- ifelse(test, "test", NA) # append to outs_root for output dir. if NULL, code uses YYYY_MM_DD

## options for simulation. the factorial selection of all combos
## comprises the set of experiments, each which are run n.rep times

num.obs <- c(16, 36, 49, 100, 400) # mean number of obs per gp
grp.var <- c(.1, 0.5, 1) # var of group random effect
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
experiments <- setDT(expand.grid(num.obs, grp.int, grp.var, n.rep))
colnames(experiments) <- c('num.obs', 'grp.int', 'grp.var', 'n.rep')
experiments[, exp := 1:.N]
fwrite(experiments, file = file.path(outs_dir, "experiments.csv"))

# make a long data.table to store the truth and the estimates
n.exp <- nrow(experiments)
n.grp <- 37 # same as nigeria
n.par <- n.grp + 2 # total effects + intercept, grp var,
n.mth <- 2 # number of methods: inla and tmb
total.rows <- n.exp * n.rep * n.par * n.mth

# nominal coverages we will check:
cp <- c(25, 50, 80, 90, 95)

master.res <- data.table(exp = rep(1:n.exp, each = total.rows / n.exp),
                         iter = rep(rep(1:n.rep, each = n.mth * n.par), n.exp),
                         seed = rep(-999.9, total.rows),
                         par = rep(c(glue("grp{1:n.grp}"), "int", "g.tau"),
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
  grp.int <- experiments[e, grp.int]
  grp.var <- experiments[e, grp.var]
  grp.tau <- 1 / grp.var

  # for each iter, we simulate data, fit both INLA and TMB, and then log the results
  for(i in 1:n.rep){

    for(mm in 1:20){
      message(glue("\n~~~~> exp {e} ({n.exp}) + iter {i} ({n.rep}): {(e-1)*n.rep+i} of {n.exp*n.rep} <~~~~\n"))
    }

    # (re)load the master res
    master.res <- fread(file.path(outs_dir, "master_results.csv"))
    # pull the results matrix for this run
    local.res <- master.res[exp == e & iter == i, ]
    rm(master.res)

    ## 0, we generate a seed by hashing a string unique to this experiment and iteration

    e.i <- glue("e{e}i{i}")

    hashed.seed <- strtoi(e.i, 36) %% (2 ^ 31 - 1) # modulo max seed integer

    local.res[exp == e & iter == i, seed := hashed.seed]
    set.seed(hashed.seed)

    ######################
    ## 1, simulate data ##
    ######################

    grp.re <- rnorm(n = n.grp, mean = 0, sd = sqrt(grp.var))
    grp.means <- grp.re + grp.int

    # true rate
    t.rate <- exp(grp.means)

    # sampled pop for each region
    s.pop <- rpois(n = n.grp, lambda = num.obs)

    # data obs
    y.obs <- rpois(n = n.grp, lambda = t.rate * s.pop)

    # save the sim obs and field
    # the obs and truth are the same for inla and tmb, so rep(..., each=2)
    local.res[grep("grp", par), truth := rep(grp.means, 2)]
    local.res[grep("int", par), truth := grp.int]
    local.res[grep("g.tau", par), truth := grp.tau]
    local.res[grep("grp", par), pop := rep(s.pop, each = 2)]
    local.res[grep("grp", par), obs := rep(y.obs, each = 2)]

    #################
    ## 2, fit INLA ##
    #################

    ## prep for INLA

    # make fixed effects design matrix
    i.data <- list(y = y.obs, E = s.pop, grp = 1:n.grp)

    # set prior for the intercept
    # alpha ~ N(mean, prec)
    priors  <- list(
      alpha = list(prior = 'normal', params = c(0, 1 / 25)) # mean, prec
    )

    # Set prior on $\text{\tau}$ s.t. the prior no scaling stdev is
    # sigma ~ N(0, prec) with sigma > 0 (truncated)

    hyper_priors <-  list(
      prec = list(prior = "logtgaussian", params=c(0, 1/25)) # mean, prec
    )

    # BYM2 model, no intercept, no covs
    i.formula <- y ~ 1 + f(grp, model = "iid",
                           hyper = hyper_priors)

    ## fit INLA
    i.fit <- hush(inla(i.formula,
                       family = "poisson", E = E,
                       data = i.data,
                       control.fixed = list(mean.intercept = 0, prec.intercept = 1/25),
                       control.predictor = list(compute=FALSE),
                       control.inla = list(strategy="simplified.laplace", fast=FALSE,
                                           int.strategy = "ccd"),
                       control.compute = list(config = TRUE), # allow for bias correct
                       verbose=FALSE, debug=TRUE, silent=TRUE))

    ################
    ## 3, fit TMB ##
    ################

    # prep for TMB
    # TMB::compile("re_iid.cpp")

    # load the model
    dyn.load( dynlib("re_iid") )

    # the data
    t.data <- list(N = n.grp,
                   y_i = y.obs,
                   e_i = s.pop,
                   options = list(adreport_on = 1,
                                  normal_trick = 0,
                                  verbose = 1),
                   flag = 1, ## only used if options[['normal_trick']]==1
                   alpha_pri = priors[['alpha']][['params']],
                   g_tau_pri = hyper_priors[['prec']][['params']])

    # starting vals
    t.params <- list(alpha = 0,
                     log_g_tau = 1,
                     g_re_i = rep(0, n.grp))

    # rand effs
    t.rand <- c('g_re_i')

    ## and lastly, we can NULL out params that are in the c++ template but
    ## which won't be used in this run. this allows us to build up a more
    ## complicated template that can take different options. named items
    ## in the list which as set to factor(NA) will be left out
    ADmap <- list()

    ## fit in TMB

    # make AD functional
    obj <- hush(MakeADFun(data       = t.data,
                          parameters = t.params,
                          random     = t.rand,
                          map        = ADmap,
                          hessian    = TRUE,
                          DLL        = "re_iid"))

    # call the optimizer
    opt0 <- hush(try(do.call("nlminb",
                             list(start     = obj$par,
                                  objective = obj$fn,
                                  gradient  = obj$gr,
                                  lower     = rep(-20, length(obj$par)),
                                  upper     = rep( 20, length(obj$par)),
                                  control   = list(trace=1)))))

    # report the estimates
    SD0 <- hush(TMB::sdreport(obj,
                              getJointPrecision = TRUE,
                              getReportCovariance = TRUE,
                              bias.correct = TRUE,
                              bias.correct.control = list(sd = TRUE)))

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
    i.draws <- hush(inla.posterior.sample(n.draws, i.fit, use.improved.mean = FALSE))

    # estimates of group totals (alpha + grp re)
    i.grp.idx <- tail(grep("Predictor", rownames(i.draws[[1]]$latent)), n.grp)
    e.i.grp   <- data.table(do.call("cbind",
                                    lapply(i.draws,
                                           function(x){x$latent[i.grp.idx]})) )
    colnames(e.i.grp) <- glue("draw{1:n.draws}")
    e.i.grp[, par := as.character(glue("grp{1:n.grp}"))][, method := "inla"]

    # inla param estimates: alpha, grp prec
    i.param.idx <- grep("Intercept", rownames(i.draws[[1]]$latent))
    i.hyper.idx <- c(grep("Precision", names(i.draws[[1]]$hyperpar)))
    e.i.param   <- data.table(do.call("cbind",
                                      lapply(i.draws,
                                             function(x){c(x$latent[i.param.idx],
                                                           x$hyperpar[i.hyper.idx])})))
    colnames(e.i.param) <- glue("draw{1:n.draws}")
    e.i.param[, par := c("int", "tau")][, method := "inla"]

    # all estimates from INLA
    e.i <- rbind(e.i.grp, e.i.param)

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # organize TMB estimates ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~

    ## take TMB draws
    mu <- c(SD0$par.fixed, SD0$par.random)

    # first, try to get the cholesky decomp of the precision
    L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)

    #if(class(L) != "try-error"){

    # take the draws
    t.draws <- rmvnorm_prec(mu = mu,
                            prec = SD0$jointPrecision,
                            n.sims = n.draws,
                            sumtozero = FALSE)

    ## separate out the tmb_draws
    t.parnames <- c(names(SD0$par.fixed), names(SD0$par.random))

    t.param.idx <- grep("alpha", t.parnames)
    t.hyper.idx <- c(grep("log_g_tau", t.parnames))
    e.t.param   <- as.matrix(t.draws[c(t.param.idx, t.hyper.idx), ])
    # transform to phi and tau
    e.t.param[2, ] <- exp(e.t.param[2, ])
    e.t.param <- data.table(e.t.param)

    colnames(e.t.param) <- glue("draw{1:n.draws}")
    e.t.param[, par := c("int", "tau")][, method := "tmb"]

    # estimates of tmb grp.res
    t.grp.idx <- grep("g_re_i", t.parnames)
    e.t.grp   <- data.table(as.matrix(t.draws[t.grp.idx, ] ))

    # add on intercept by draw
    e.t.grp <- data.table( t(t(e.t.grp) + t.draws[t.param.idx, ]) )

    colnames(e.t.grp) <- glue("draw{1:n.draws}")
    e.t.grp[, par := as.character(glue("grp{1:n.grp}"))][, method := "tmb"]

    # all estimates from TMB
    e.t <- rbind(e.t.grp, e.t.param)

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # summarize INLA and TMB ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~

    est <- rbind(e.i, e.t)

    # get est := mean(draws)
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
