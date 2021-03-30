# this script launches discrete simulations to compare fits from TMB and R-INLA
# aoz - 2020

###############
# USER INPUTS #
###############

test <- FALSE

## root dirs

data_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects'
code_root <- '/home/merms/Documents/GitRepos/tmb_inla_comp/'
outs_root <- '/home/merms/Documents/Research/2020_tmb_v_inla/discrete_sim'
outs_name <- ifelse(test, "test", NA) # append to outs_root for output dir. if NULL, code uses YYYY_MM_DD

## options for simulation. the factorial selection of all combos
## comprises the set of experiments, each which are run n.rep times

num.obs <- c(16, 36, 49, 100, 400) # mean number of obs per area
bym.var <- 0.5 # total var of BYM2 effect
bym.phi <- c(0.25, .5, .75, .9) # mixing param determining amount of BYM which is structured
bym.int <- -3 # intercept of bym2 GMRF
n.rep   <- 25 # number of reps per experiments

# set number of draws for summarization of fits
n.draws <- 500

#############
# SETUP ENV #
#############

setwd(file.path(code_root, "discrete_sim"))

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
source(file.path(code_root, 'discrete_sim/discrete_sim_funcs.R'))

# load prepped nigeria data
load(file.path(data_root, "prepped_discrete_obj.RData"))

# convert to sf obj
nga1 <- st_as_sf(nga1)

# grab the adjacency matrix for NGA
adj.mat <- W.nga1

# prep the scaled ICAR precision we need for the BYM2 model
Q.icar <- make_BYM2_struct_scaled_prec(adj.mat,
                                       sparse.mat = TRUE)


#######################################################
# make master objects for saving truths and estimates #
#######################################################

# make the set of experiments we will run
experiments <- setDT(expand.grid(num.obs, bym.int, bym.var, bym.phi, n.rep))
colnames(experiments) <- c('num.obs', 'bym.int', 'bym.var', 'bym.phi', 'n.rep')
experiments[, exp := 1:.N]
fwrite(experiments, file = file.path(outs_dir, "experiments.csv"))

# make a long data.table to store the truth and the estimates
n.exp <- nrow(experiments)
n.reg <- nrow(nga1)
n.par <- n.reg + 3 # total effects + intercept, phi, tau
n.mth <- 2 # number of methods: inla and tmb
total.rows <- n.exp * n.rep * n.par * n.mth

# nominal coverages we will check:
cp <- c(25, 50, 80, 90, 95)

master.res <- data.table(exp = rep(1:n.exp, each = total.rows / n.exp),
                         iter = rep(rep(1:n.rep, each = n.mth * n.par), n.exp),
                         seed = rep(-999.9, total.rows),
                         par = rep(c(glue("bym{1:n.reg}"), "int", "phi", "tau"),
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
  bym.var <- experiments[e, bym.var]
  bym.tau <- 1 / bym.var
  bym.phi <- experiments[e, bym.phi]
  bym.int <- experiments[e, bym.int]

  # for each iter, we simulate data, fit both INLA and TMB, and then log the results
  for(i in 1:n.rep){

    message("\n\n\n\n\n\n\n\n\n\n")
    for(mm in 1:10){
      message(glue("\n~~~~> exp {e} ({n.exp}) + iter {i} ({n.rep}): {(e-1)*n.rep+i} of {n.exp*n.rep} <~~~~\n"))
    }
    message("\n\n\n\n\n\n\n\n\n\n")

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

    # make precision for bym2 | prec, phi
    Q.bym <- make_BYM2_joint_unstruct_scaled_struc_prec_mat(adj.mat = adj.mat,
                                                            phi = bym.phi,
                                                            tau = 1 / bym.var,
                                                            sparse.mat = TRUE)

    # sim true bym2 field
    bym <- rbym2_simul_1constr_prec(mu = 0, # intercept is passed into sim_field_data()
                                    prec = Q.bym,
                                    n.sims = 1,
                                    sumtozero = TRUE)$x.c # take the constrained draw
    bym.total <- bym[1:n.reg]
    bym.struc <- bym[(n.reg + 1):(2 * n.reg)]

    # sampled pop for each region
    s.pop <- rpois(n = n.reg, lambda = num.obs)

    # sim data
    obs.risk <- sim_field_data(pop = s.pop,
                               int = bym.int,
                               field = bym.total,
                               lik = "poisson")

    # data obs
    y.obs <- obs.risk[, obs]

    # true risk
    t.risk <- obs.risk[, risk]

    # save the sim obs and field
    # the obs and truth are the same for inla and tmb, so rep(..., each=2)
    local.res[grep("bym", par), truth := rep(bym.total, 2)]
    local.res[grep("int", par), truth := bym.int]
    local.res[grep("phi", par), truth := bym.phi]
    local.res[grep("tau", par), truth := bym.tau]
    local.res[grep("bym", par), pop := rep(s.pop, each = 2)]
    local.res[grep("bym", par), obs := rep(y.obs, each = 2)]

    #################
    ## 2, fit INLA ##
    #################

    ## prep for INLA

    i.data <- list(
      y = y.obs, E = s.pop, region = nga1$ADM1_CODE
    )

    # Get the adjacency matrix INLA needs (1s on diagonal)
    i.adj.mat <- adj.mat + Diagonal(nrow(adj.mat))

    # set prior for the intercept
    # alpha ~ N(mean, prec)
    priors  <- list(
      alpha = list(prior = 'normal', params = c(0, 1 / 25)) # mean, prec
    )

    # Set prior on $\text{logit}(\phi)$ s.t. the mixing param prior is
    # phi ~ beta(1/2,1/2)$

    # Set prior on $\text{\tau}$ s.t. the prior no scaling stdev is
    # sigma ~ N(0, prec) with sigma > 0 (truncated)

    hyper_priors <-  list(
      phi  = list(prior = "logitbeta", params=c(0.5, 0.5)),
      prec = list(prior = "logtgaussian", params=c(0, 1/25)) # mean, prec
    )

    # BYM2 model, no intercept, no covs
    i.formula <- y ~ 1 + f(region, model = "bym2",
                           graph = i.adj.mat,
                           hyper = hyper_priors,
                           constr = TRUE)

    ## fit INLA
    i.fit <- inla(i.formula,
                  family = "poisson", E = E,
                  data = i.data,
                  control.fixed = list(mean.intercept = 0, prec.intercept = 1/25),
                  control.predictor = list(compute=FALSE),
                  control.inla = list(strategy="simplified.laplace", fast=FALSE,
                                      int.strategy = "ccd"),
                  control.compute = list(config = TRUE), # allow for bias correct
                  verbose=FALSE, debug=TRUE, silent=TRUE)

    ################
    ## 3, fit TMB ##
    ################

    ## prep for TMB

    # load the model
    dyn.load( dynlib("bym2") )

    # the data
    t.data <- list(N = n.reg,
                   y_i = y.obs,
                   e_i = s.pop,
                   Q_struc = Q.icar,
                   options = list(adreport_on = 1,
                                  normal_trick = 0,
                                  verbose = 0),
                   flag = 1, ## only used if options[['normal_trick']]==1
                   bym2_alpha_pri = priors[['alpha']][['params']],
                   bym2_phi_pri = hyper_priors[['phi']][['params']],
                   bym2_tau_pri = hyper_priors[['prec']][['params']])

    # starting vals
    t.params <- list(alpha = 0,
                     logit_phi = 0,
                     log_tau = 0,
                     Epsilon_s = matrix(0, ncol = 1, nrow = 2 * nrow(t.data[['Q_struc']])))

    # rand effs
    t.rand <- c('Epsilon_s')

    ## and lastly, we can NULL out params that are in the c++ template but
    ## which won't be used in this run. this allows us to build up a more
    ## complicated template that can take different options. named items
    ## in the list which as set to factor(NA) will be left out
    ADmap <- list()

    ## fit in TMB

    # make AD functional
    obj <- MakeADFun(data       = t.data,
                     parameters = t.params,
                     random     = t.rand,
                     map        = ADmap,
                     hessian    = TRUE,
                     DLL        = "bym2")

    # call the optimizer
    opt0 <- try(do.call("nlminb",
                        list(start     = obj$par,
                             objective = obj$fn,
                             gradient  = obj$gr,
                             lower     = rep(-10, length(obj$par)),
                             upper     = rep( 10, length(obj$par)),
                             control   = list(trace=1))))

    # report the estimates
    SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                         getReportCovariance = TRUE,
                         bias.correct = TRUE,
                         bias.correct.control = list(sd = TRUE))

    ## summary(SD0, report = 'report')

    R0 <- obj$report()

    ##########################
    ## 4, Summarize results ##
    ##########################

    # we want:
    #  the field estimates and their coverage of the truth
    #  the param estimates and their coverage of the truth

    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    # organize INLA estimates ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~

    # take INLA draws
    i.draws <- inla.posterior.sample(n.draws, i.fit, use.improved.mean = FALSE)

    #estimates of structured component
    i.struc.idx <- tail(grep("region", rownames(i.draws[[1]]$latent)), n.reg)
    e.i.struc   <- data.table(do.call("cbind",
                                      lapply(i.draws,
                                             function(x){x$latent[i.struc.idx]})) )
    colnames(e.i.struc) <- glue("draw{1:n.draws}")
    e.i.struc[, par := as.character(glue("bym{1:n.reg}"))][, method := "inla"]


    # estimates of inla field
    # i.field.idx <- grep("Predictor", rownames(i.draws[[1]]$latent))
    i.total.idx <- grep("region", rownames(i.draws[[1]]$latent))[1:n.reg]
    e.i.total   <- data.table(do.call("cbind",
                                      lapply(i.draws,
                                             function(x){x$latent[i.total.idx]})) )
    colnames(e.i.total) <- glue("draw{1:n.draws}")
    e.i.total[, par := as.character(glue("bym{1:n.reg}"))][, method := "inla"]

    # inla param estimates: alpha, phi, prec
    i.param.idx <- grep("Intercept", rownames(i.draws[[1]]$latent))
    i.hyper.idx <- c(grep("Phi", names(i.draws[[1]]$hyperpar)),
                     grep("Precision", names(i.draws[[1]]$hyperpar)))
    e.i.param   <- data.table(do.call("cbind",
                                      lapply(i.draws,
                                             function(x){c(x$latent[i.param.idx],
                                                           x$hyperpar[i.hyper.idx])})))
    colnames(e.i.param) <- glue("draw{1:n.draws}")
    e.i.param[, par := c("int", "phi", "tau")][, method := "inla"]

    # all estimates from INLA
    e.i <- rbind(e.i.total, e.i.param)

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # organize TMB estimates ~
    #~~~~~~~~~~~~~~~~~~~~~~~~~

    ## take TMB draws

    mu <- c(SD0$par.fixed, SD0$par.random)

    # first, try to get the cholesky decomp of the precision
    L <- try(suppressWarnings(Cholesky(SD0$jointPrecision, super = T)), silent = TRUE)

    #if(class(L) != "try-error"){

    # take the draws
    t.draws <- rbym2_simul_1constr_prec(mu = mu,
                                        prec = SD0$jointPrecision,
                                        n.sims = n.draws,
                                        sumtozero = TRUE,
                                        constrain.idx = tail(1:length(mu), n.reg)) # constrain the last n.reg

    #}

    # take the constrained draws
    t.draws <- t.draws$x.c

    ## separate out the tmb_draws
    t.parnames <- c(names(SD0$par.fixed), names(SD0$par.random))

    t.param.idx <- grep("alpha", t.parnames)
    t.hyper.idx <- c(grep("phi", t.parnames),
                     grep("tau", t.parnames))
    e.t.param   <- as.matrix(t.draws[c(t.param.idx, t.hyper.idx), ])
    # transform to phi and tau
    e.t.param[2, ] <- plogis(e.t.param[2, ])
    e.t.param[3, ] <- exp(e.t.param[3, ])
    e.t.param <- data.table(e.t.param)

    colnames(e.t.param) <- glue("draw{1:n.draws}")
    e.t.param[, par := c("int", "phi", "tau")][, method := "tmb"]

    # estimates of tmb field
    t.total.idx <- head(grep("Epsilon", t.parnames), n.reg)
    e.t.total   <- data.table(as.matrix(t.draws[t.total.idx, ] ))

    # add on intercept by draw
    # e.t.field <- data.table( t(t(e.t.total) + t.draws[t.param.idx, ]) )

    colnames(e.t.total) <- glue("draw{1:n.draws}")
    e.t.total[, par := as.character(glue("bym{1:n.reg}"))][, method := "tmb"]

    # all estimates from TMB
    e.t <- rbind(e.t.total, e.t.param)

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
