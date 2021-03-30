## ##########################################
## setup the environment for singularity R ##
## ##########################################

message('---- ON SCRIPT 2: setup experiment settings')

## to help with loading in arguments in loopvars
options(stringsAsFactors = FALSE)

## Set core_repo location and tmb_repo loc
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)

## grab libraries and functions from MBG code
setwd(core_repo)
commondir    <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

## Load MBG packages and functions
source(paste0(core_repo, 'mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(scales)
library(RandomFields)

## set the location of the pardiso solver
inla.setOption("pardiso.license", "~/sys/licenses/pardiso.lic")

## Now we can switch to the TMB repo
setwd(tmb_repo)
source('./realistic_sim_utils.R')

## set other needed global args 
modeling_shapefile_version <- '2019_10_10'

## setup is now done. setup some parameters for this simulation

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################################
## load in the loopvars from launch and setup all the params for this job ##
############################################################################

# these three lines were moved to script 1 when switching to array jobs on 08mar2020
# main.dir <- sprintf('/ihme/scratch/users/azimmer/tmb_inla_sim/%s', main.dir.name)
# out.dir  <- sprintf('%s/%06d', main.dir, exp.lvid)
# loopvars <- read.csv(file = paste0(main.dir, '/loopvars.csv'))

## create a directory for some common objects that can be shared by all experiments launched
common.dir <- sprintf('%s/common/', main.dir)
dir.create(common.dir, recursive = TRUE, showWarnings = F)

## create a write a file to keep track of job status
jobtrack.dir <- paste0(common.dir, 'job_tracking/')
dir.create(jobtrack.dir, recursive = TRUE, showWarnings = F)
write.table(x=matrix(c(sim.loop.ct, 2), ncol=2), 
            file = paste0(jobtrack.dir, 
                          sprintf('exp_%06d_iter_%06d.csv', exp.lvid, exp.iter)), sep=',', 
            row.names=F)

## create some directories for output organization
dir.create(sprintf('%s/simulated_obj', out.dir), recursive = TRUE, showWarnings = F)
dir.create(sprintf('%s/modeling/inputs', out.dir), recursive = TRUE, showWarnings = F)
dir.create(sprintf('%s/modeling/outputs/tmb', out.dir), recursive = TRUE, showWarnings = F)
dir.create(sprintf('%s/modeling/outputs/inla', out.dir), recursive = TRUE, showWarnings = F)
dir.create(sprintf('%s/validation', out.dir), showWarnings = F)

## load in all parameters for this experiment
reg             <- as.character(loopvars[exp.lvid, 1])
year_list       <- eval(parse(text = loopvars[exp.lvid, 2]))
cov_names       <- eval(parse(text = as.character(loopvars[exp.lvid, 3])))
cov_measures    <- eval(parse(text = as.character(loopvars[exp.lvid, 4])))
betas           <- eval(parse(text = as.character(loopvars[exp.lvid, 5])));if(mean(is.na(betas))==1) betas <- NULL

alpha           <- as.numeric(loopvars[exp.lvid, 6]);if(is.na(alpha)) alpha <- NULL
sp.range        <- as.numeric(loopvars[exp.lvid, 7])
sp.var          <- as.numeric(loopvars[exp.lvid, 8])
sp.alpha        <- as.numeric(loopvars[exp.lvid, 9])
clust.var       <- as.numeric(loopvars[exp.lvid, 10]);if(is.na(clust.var)) clust.var <- NULL

t.rho           <- as.numeric(loopvars[exp.lvid, 11])
mesh_s_params   <- as.character(loopvars[exp.lvid, 12])
n.clust         <- as.numeric(loopvars[exp.lvid, 13])
m.clust         <- as.numeric(loopvars[exp.lvid, 14])
sample.strat    <- eval(parse(text = as.character(loopvars[exp.lvid, 15])))
obs.loc.strat   <- sample.strat[['obs.loc.strat']]
urban.pop.pct   <- sample.strat[['urban.pop.pct']]
urban.strat.pct <- sample.strat[['urban.strat.pct']]

cores           <- as.numeric(loopvars[exp.lvid, 16]) 
ndraws          <- as.numeric(loopvars[exp.lvid, 17])
alphaj.pri      <- eval(parse(text = as.character(loopvars[exp.lvid, 18]))) ## normal mean and sd ## TODO pass this to INLA and TMB
clust.prec.pri    <- eval(parse(text = as.character(loopvars[exp.lvid, 19])))  ## gamma for clust preciion with shape and inv-scale ## TODO pass this to INLA and TMB
inla.int.strat  <- as.character(loopvars[exp.lvid, 20]) ## can be one of: 'eb', 'ccd', 'grid'

inla.approx     <- as.character(loopvars[exp.lvid, 21]) ## can be 'gaussian', 'simplified.laplace' (default) or 'laplace'
l.tau.pri     <- NULL  ## taken from INLA spde mesh obj
l.kap.pri     <- NULL  ## taken from INLA spde mesh obj
Nsim          <- as.numeric(loopvars[exp.lvid, 22]) ## number of times to repeat simulation
data.lik      <- as.character(loopvars[exp.lvid, 23])
norm.var      <- as.numeric(loopvars[exp.lvid, 24])
norm.prec.pri <- eval(parse(text = as.character(loopvars[exp.lvid, 25])))

bias.correct <- as.logical(loopvars[exp.lvid, 26])
sd.correct   <- as.logical(loopvars[exp.lvid, 27])
matern.pri   <- eval(parse(text = as.character(loopvars[exp.lvid, 28]))) ## pc prior for matern
fix.locs     <- as.logical(loopvars[exp.lvid, 29]) ## fix locations across experiments
fix.gp       <- as.logical(loopvars[exp.lvid, 30]) ## fix locations across experiments


## TODO? add in some validation options? or maybe just always do them all

## make a vector of all possible params to use throughout (mostly in validation)
true.params <- data.table(param = c('int',
                                    cov_names,
                                    'sp gp range',
                                    'sp gp var',
                                    'clust var',
                                    'norm obs var'
                                    ),
                          truth = c(ifelse(is.null(alpha), NA, alpha),
                                    switch(2 - is.null(betas), rep(NA, length(cov_names)), betas),
                                    sp.range,
                                    sp.var,
                                    ifelse(is.null(clust.var), NA, clust.var),
                                    ifelse(data.lik == 'normal', norm.var, NA)
                                    )
                          )

write.table(file = sprintf('%s/simulated_obj/true_param_table.csv', out.dir),
            x = true.params, sep=',' ,row.names = FALSE)

## from these imputs, make a table of covariate names and measures
covs <- data.table(name = cov_names, meas = cov_measures)

## I hardcode a few other options that are useful sometimes when running interactively
## these can probably be deleted...

## to make it easier to run real data from this code, usually have SIM=TRUE, REAL=FALSE
use_real_data <- FALSE
use_sim_data  <- TRUE

## should we save this data to mbg input_data? useful if we want to run sim data through lbd_core
save.as.input.data <- FALSE
data.tag <- '_allyrs_clust'

## end of user inputs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## transform some inputs into other useful quantities

## convert sp.range and sp.var into sp.kappa for rspde function
sp.kappa   <- sqrt(8) / sp.range
logkappa   <- log(sp.kappa)
sp.tau     <- sqrt(1 / (4 * pi * sp.kappa ^ 2 * sp.var)) ## sp.var = 1/(4*pi*kappa^2*tau^2)
logtau     <- log(sp.tau)
trho_trans <- log((-1 - t.rho) / (t.rho - 1))

###########################################
## load in region/counry shapes and covs ##
###########################################

## load in the region shapefile and prep the boundary
if(!file.exists(sprintf('%s/poly_shape.rdata', common.dir))){
  
  gaul_list           <- get_adm0_codes(reg)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]
  rm(simple_polygon_list)

  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]
  rm(raster_list)

  ## save these since they will be common across runs (assuming the geography is fixed)
  if(exp.lvid == 1){
    ## only save from the first one to avoid simultaneously writing to the same file
    save(gaul_list,
         subset_shape, 
         simple_polygon,
         simple_raster,
         pop_raster,
         file = sprintf('%s/poly_shape.rdata', common.dir))
  }
  
}else{
  message('------ loading in pre-made simple raster, polygon and population objects')
  load(sprintf('%s/poly_shape.rdata', common.dir))
}

