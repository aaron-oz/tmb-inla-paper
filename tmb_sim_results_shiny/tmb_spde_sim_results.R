# source('/home/merms/Documents/GitRepos/tmb_inla_comp/tmb_sim_results_shiny/tmb_spde_sim_results.R')
library(shiny)
library(glue)
library(TMB)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(plyr)
library(ggplot2)
library(data.table)

# prep some and load objects
user      <- Sys.info()['user']
core_repo <- sprintf('/share/code/geospatial/%s/lbd_core/', user)
tmb_repo  <- sprintf('/homes/%s/tmb_inla_comp', user)
tmb_repo  <- "~/Documents/GitRepos/tmb_inla_comp"

setwd(tmb_repo)
source('./realistic_sim_utils.R')

main.dir.name <- "spatial-stats-sub-2020_09_05_21_56_12" # using RandomFields to sim GP
main.dir.root <- "/ihme/scratch/users/azimmer/tmb_inla_sim"
main.dir.root <- "/home/merms/Documents/Research/2020_tmb_v_inla/tmb_inla_sim"

main.dir <- file.path(main.dir.root, main.dir.name)
compar.dir <- file.path(main.dir, "comparisons")

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


## reload the prepped file
message('--reloading in combined summary metrics')
summary.metrics <- fread(sprintf('%sall_summary_metrics.csv', compar.dir))

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

## convert mesh params params to low/med/high res
summary.metrics[st_mesh_nodes == 3631, mesh.res := 'low']
summary.metrics[st_mesh_nodes == 7922, mesh.res := 'med']
summary.metrics[st_mesh_nodes == 13869, mesh.res := 'high']
## set the order
summary.metrics$mesh.res <- factor(summary.metrics$mesh.res, levels = c('low', 'med', 'high'))

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## bias plots prep
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                                       'fit_type', 'clust.var.cat', 'mesh.res',
                                       "inla.int.strat", "inla.approx"), with=F]
fe.mean.long <- melt(key.summ.metrics, id=c('data.lik', 'norm.var', 'n.clust',
                                            'fit_type', 'clust.var.cat', 'mesh.res',
                                            "inla.int.strat", "inla.approx"))

## drop NAs. eg there is no gauss_prec_bias in a binom model,
## and there is no clust_prec_bias when clust.var==NA
fe.mean.long <- fe.mean.long[!is.na(value),]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## write the overview text
str1 <- "<b>Overview</b>"
str2 <- "These plots show results from the continuous spatial Gaussian process (GP) simulation experiments describedby A. Osgood-Zimmerman and J. Wakefield in <i>A Statistical Introduction to Template Model Builder: A Flexible Tool for Spatial Modeling</i><sup>1</sup>."
str3 <- 'Briefly, continuous Gaussian fields with Matern covariance kernels were simulated over the spatial domain of Nigeria. The Stochastic Process Differential Equation (SPDE) representation<sup>2</sup> is used to approximate the GPs during inference using both Template Model Builder (TMB<sup>3</sup>) and Integrated Nested Laplace Approximations (INLA<sup>4</sup>), and comparisons between the approximate posterior distributions and the truth are made for 16,128 experimental levels each replicated 25 times.'
str4 <- 'These plots let you explore the results while varying<ul>
      <li>the data likelihood (Gaussian or binomial),</li>
      <li>the resolution of the SPDE mesh (low=3631, med=7922, high=13869 vertices, respectively),</li>
      <li>the INLA approximations (Full (L)aplace, (S)implified Laplace, and (G)aussian) and hyperparameter integration strategies ( (C)CD or (E)mpirical Bayes),</li>
      <li>and the quantile range of the Monte Carlo error bars.</li>
    </ul>
Adjust the height of the plots until they are comfortable to view on your screen'
overview <- HTML(paste(str1, str2, str3, str4, sep = '<br/>'))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define app pages
ui <- fluidPage(

  # App title ----
  titlePanel("TMB and R-INLA SPDE Simulation Results"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for input ----
    sidebarPanel(

      # Input: checkboxes for methods to compare
      checkboxGroupInput("methodType",
                         #                         h3("Methods"),
                         label = "Select all methods to compare",
                         choices = list("TMB" = 1,
                                        "INLA + G + EB" = 2,
                                        "INLA + G + CCD" = 3,
                                        "INLA + SL + EB" = 4,
                                        "INLA + SL + CCD" = 5,
                                        "INLA + L + EB" = 6,
                                        "INLA + L + CCD" = 7
                                        ),
                         selected = c(1, 7) # initial selected options
                         ),

      # Input: toggle binomial or gaussian data
      radioButtons("dataType",
                   label = "Select data likelihood",
                   choices = list("Binomial" = 1,
                                  "Gaussian" = 2),
                   selected = c(1) # initial selected options
                   ),

      # Input: select Gaussian data observation variance
      # TODO: grey out if binomial data type is selected
      ## radioButtons("obsVar",
      ##              label = "Select Gaussian observation variance (only needed for some plots)",
      ##              choices = list("0.1^2" = 1,
      ##                             "0.2^2" = 2,
      ##                             "0.4^2" = 3),
      ##              selected = c(1) # initial selected options
      ##              ),

      ## Input: select Gaussian data observation variance
      radioButtons("meshRes",
                   label = "Select SPDE triangulation resolution",
                   choices = list("low" = 1,
                                  "med"  = 2,
                                  "high" = 3,
                                  "all"  = 4),
                   selected = c(4) # initial selected options
                   ),

      helpText("Note: when the 'all' mesh resolution is picked,",
               "the results will be averaged over the 'low',",
               "'med', and 'high' resolutions."),

      # Input: Slider for the monte carlo error bar width/quantiles
      sliderInput(inputId = "coverage",
                  label = "Monte Carlo Quantile Coverage (equal-tailed intervals)",
                  min = 25,
                  max = 100,
                  step = 5,
                  value = 50),

      # Input: Slider for the plot height
      sliderInput(inputId = "plotHeight",
                  label = "Adjust plot height (pixels)",
                  min = 400,
                  max = 2000,
                  step = 100,
                  value = 800),

      helpText(HTML("<sup>1</sup> <i>Forthcoming</i>")),

      helpText(HTML("<sup>2</sup> Lindgren, Finn, Johan Lindström, and Håvard Rue. <i>An explicit link between Gaussian fields and Gaussian Markov random fields: The SPDE approach</i>. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73: 423-498. https://doi.org/10.1111/j.1467-9868.2011.00777.x")),

      helpText(HTML("<sup>3</sup> Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, & Bradley M. Bell. <i>TMB: Automatic Differentiation and Laplace Approximation.</i> Journal of Statistical Software [Online], 70.5 (2016): 1 - 21. Web. 16 Mar. 2021
")),

helpText(HTML("<sup>4</sup> Rue, Håvard, Sara Martino, and Nicolas Chopin. <i>Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations.</i> Journal of the Royal Statistical Society: Series B (Statistical Methodology) 71.2 (2009): 319-392."))

),

# Main panel for displaying outputs ----
mainPanel(
  tabsetPanel(
    tabPanel("Parameter Bias",
             div(overview, style = "border:3px;
                                        border-style:solid;
                                        border-radius:8px;
                                        border-color:#689CBE;
                                        padding: 1em;"),
             htmlOutput("biasPlotText", style = "border:3px;
                                                     border-style:solid;
                                                     border-radius:8px;
                                                     border-color:Gray;
                                                     padding: 1em;"),
             plotOutput("biasPlot")
             ),
    tabPanel("Total Field Coverage",
             div(overview, style = "border:3px;
                                        border-style:solid;
                                        border-radius:8px;
                                        border-color:#689CBE;
                                        padding: 1em;"),
             htmlOutput("fieldPlotText", style = "border:3px;
                                                      border-style:solid;
                                                      border-radius:8px;
                                                      border-color:Gray;
                                                      padding: 1em;"),
             plotOutput("fieldPlot"),
             )
    #        tabPanel("Decile Field Coverage", plotOutput("decilePlot"))
  )
)
  )
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define server logic required to make plots
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when input (input$coverage) change
  # 2. Its output type is a plot
  output$fieldPlot <- renderPlot({

    cov.width       <- input$coverage / 100
    dl              <- c("binom", "normal")[as.numeric(input$dataType)]
    model.selection <- as.numeric(input$methodType)
    mesh.resolution <- as.numeric(input$meshRes)
    ## norm.var        <- (c(0.1, 0.2, 0.4) ^ 2)[as.numeric(input$obsVar)]

    models.to.plot <- c("TMB",
                        "INLA_EB_G",
                        "INLA_CCD_G",
                        "INLA_EB_S",
                        "INLA_CCD_S",
                        "INLA_EB_L",
                        "INLA_CCD_L")[model.selection]

    model.plot.names <- c("TMB",
                          "INLA: G + EB",
                          "INLA: G + CCD",
                          "INLA: SL + EB",
                          "INLA: SL + CCD",
                          "INLA: L + EB",
                          "INLA: L + CCD")[model.selection]


    ## make the plot title
    pt <- sprintf('Average Pixel Coverage across Spatial Domain for %s Observations',
                  stringr::str_to_title(dl))
    #    if(dl == "normal"){
    #      pt <- paste0(pt, sprintf(" with %0.02f Variance:\n", norm.var))
    #    }else{
    pt <- paste0(pt, ":\n")
    #    }
    pt <- paste(pt, sprintf('Median and %0.0f%% Quantile Bars across Monte Carlo Replicates', cov.width*100))

    # subset to all models desired
    long.cov.sub <- rbindlist(lapply(models.to.plot, function(x){subset(long.cov,
                                                                        data.lik == dl &
                                                                          grepl(x, long.cov$fit_type))}))

    ## subset to mesh
    if(mesh.resolution %in% c('low', 'med', 'high')){
      long.cov.sub <- subset(long.cov.sub, mesh.res == mesh.resolution)
    } ## ow, it is 'all', and we average across it

    # nice plot legend names
    for(nn in 1:length(model.plot.names)){
      long.cov.sub[fit_type == models.to.plot[nn], fit_type := model.plot.names[nn]]
    }

    ## facet by observations
    long.cov.sum <-  summarySE(long.cov.sub, measurevar="obs_cov",
                               groupvars = c('n_target_cov',
                                             'fit_type', 'n.clust',
                                             'data.lik',
                                             'clust.var.cat',
                                             'norm.var'),
                               conf.interval = cov.width)
    ## groups to plot lines
    long.cov.sum$line_group <- paste(long.cov.sum$fit_type, long.cov.sum$clust.var, sep='_')

    ## set facet names using the labeller
    ## annoying format. must be as_labeller(c(`facet_name1`='Preferred Facet Name1', ...))
    ## facet_labs <- as_labeller(eval(parse(text=paste0('c(', paste0('`', sort(unique(long.cov$n.clust)), '`', '=', '\'', paste0('Num. Cluster Obs: ', sort(unique(long.cov$n.clust))), '\'', sep='', collapse=','), ')'))))

    ## final rename for plotting
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
      theme(text = element_text(size = 20))  + # increase base size of all font
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

    text = paste("\n   The following is text that'll appear in a plot window.\n",
                 "       As you can see, it's in the plot window\n",
                 "       One might imagine useful information here")
    fit_coverage_CI_summary_text <-
      ggplot() +
      annotate("text", x = 4, y = 25, size=8, label = text) +
      theme_void()

    return(fit_coverage_CI_summary)

  },

  width = "auto",
  height = reactive(input$plotHeight))

  output$fieldPlotText <- renderUI({

    cov.perc        <- paste0(input$coverage, '%')
    dl              <- c("binomial", "Gaussian")[as.numeric(input$dataType)]
    model.selection <- as.numeric(input$methodType)
    mesh.resolution <- c('low', 'medium', 'high', 'all')[as.numeric(input$meshRes)]
    ## norm.var        <- (c(0.1, 0.2, 0.4) ^ 2)[as.numeric(input$obsVar)]

    model.plot.names <- c("TMB",
                          "INLA: G + EB",
                          "INLA: G + CCD",
                          "INLA: SL + EB",
                          "INLA: SL + CCD",
                          "INLA: L + EB",
                          "INLA: L + CCD")[model.selection]
    model.plot.names.string <- model.plot.names[1]
    model.plot.names <- model.plot.names[-1]
    while(length(model.plot.names) > 1){
      model.plot.names.string <- paste(model.plot.names.string, model.plot.names[1], sep = ', ')
      model.plot.names <- model.plot.names[-1]
    }
    if(length(model.plot.names) == 1){
      if(length(model.selection) > 2){
        model.plot.names.string <- paste0(model.plot.names.string, ',')
      }
      model.plot.names.string <- paste(model.plot.names.string, model.plot.names[1], sep = ' and ')
    }

    title <- '<b>Figure Description</b>'

    str1 <- paste0('Comparison of the average estimated field coverage of the ',
                   'simulated truth, faceted by ',
                   ifelse(dl == 'Gaussian', 'observation variances and ', ''),
                   'the number of clusters, from ',
                   model.plot.names.string,
                   ', plotted against the target nominal coverage, for ', dl,
                   ' observation experiments with ',
                   ifelse(mesh.resolution == 'all', 'low, medium, and high resolution SPDE triangulations.',
                          paste0('the ', mesh.resolution, ' resolution SPDE triangulation.')),
                   ' Colors represent different cluster (i.i.d nugget) variances',
                   ' used in an experiment. Each point',
                   ' is the median average coverage of ',
                   ifelse(mesh.resolution == 'all', 'three experiments, calculated across 75 replicates,',
                          'one experiment, calculated across 25 replicates,'),
                   '  and the bars represent the middle ', cov.perc, ' quantile range of the ',
                   'average coverage across replicates. Slight jitter has been applied to improve legibility.' )

    HTML(paste(title, str1, sep = '<br/>'))
  })

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  output$biasPlot <- renderPlot({

    cov.width       <- input$coverage / 100
    dl              <- c("binom", "normal")[as.numeric(input$dataType)]
    model.selection <- as.numeric(input$methodType)
    mesh.resolution <- as.numeric(input$meshRes)
    ## norm.var        <- (c(0.1, 0.2, 0.4) ^ 2)[as.numeric(input$obsVar)]

    models.to.plot <- c("TMB",
                        "INLA_EB_G",
                        "INLA_CCD_G",
                        "INLA_EB_S",
                        "INLA_CCD_S",
                        "INLA_EB_L",
                        "INLA_CCD_L")[model.selection]

    model.plot.names <- c("TMB",
                          "INLA: G + EB",
                          "INLA: G + CCD",
                          "INLA: SL + EB",
                          "INLA: SL + CCD",
                          "INLA: L + EB",
                          "INLA: L + CCD")[model.selection]


    ## make the plot title
    pt <- sprintf('Parameter Bias for %s Observations',
                  stringr::str_to_title(dl))
    #    if(dl == "normal"){
    #     pt <- paste0(pt, sprintf(" with %0.02f Variance:\n", norm.var))
    #  }else{
    pt <- paste0(pt, ":\n")
    # }
    pt <- paste(pt, sprintf('Median and %0.0f%% Quantile Bars across Monte Carlo Replicates', cov.width*100))

    # subset to all models desired
    long.cov.sub <- rbindlist(lapply(models.to.plot, function(x){subset(fe.mean.long,
                                                                        data.lik == dl &
                                                                          grepl(x, fe.mean.long$fit_type))}))

    ## subset to mesh
    if(mesh.resolution %in% c('low', 'med', 'high')){
      long.cov.sub <- subset(long.cov.sub, mesh.res == mesh.resolution)
    } ## ow, it is 'all', and we average across it

    # nice plot legend names
    for(nn in 1:length(model.plot.names)){
      long.cov.sub[fit_type == models.to.plot[nn], fit_type := model.plot.names[nn]]
    }

    # get intervals
    fe.mean.long.sum <-  summarySE(long.cov.sub, measurevar="value",
                                   groupvars=c("data.lik", 'norm.var', 'n.clust',
                                               "fit_type", 'clust.var.cat',
                                               'variable'),
                                   conf.interval = cov.width)

    ## final rename for plotting
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
      theme(text = element_text(size = 20))  + # increase base size of all font
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

    return(fe_bias_summary)

  },
  width = "auto",
  height = reactive(input$plotHeight))

  output$biasPlotText <- renderUI({

    cov.perc        <- paste0(input$coverage, '%')
    dl              <- c("binomial", "Gaussian")[as.numeric(input$dataType)]
    model.selection <- as.numeric(input$methodType)
    mesh.resolution <- c('low', 'medium', 'high', 'all')[as.numeric(input$meshRes)]
    ## norm.var        <- (c(0.1, 0.2, 0.4) ^ 2)[as.numeric(input$obsVar)]

    model.plot.names <- c("TMB",
                          "INLA: G + EB",
                          "INLA: G + CCD",
                          "INLA: SL + EB",
                          "INLA: SL + CCD",
                          "INLA: L + EB",
                          "INLA: L + CCD")[model.selection]
    model.plot.names.string <- model.plot.names[1]
    model.plot.names <- model.plot.names[-1]
    while(length(model.plot.names) > 1){
      model.plot.names.string <- paste(model.plot.names.string, model.plot.names[1], sep = ', ')
      model.plot.names <- model.plot.names[-1]
    }
    if(length(model.plot.names) == 1){
      if(length(model.selection) > 2){
        model.plot.names.string <- paste0(model.plot.names.string, ',')
      }
      model.plot.names.string <- paste(model.plot.names.string, model.plot.names[1], sep = ' and ')
    }

    title <- '<b>Figure Description</b>'

    str1 <- paste0('Comparison of the estimated parameter bias from ',
                   model.plot.names.string,
                   ', plotted against the number of cluster observations for the ')

    str1 <- paste0(str1, dl, ' data experiments',
                   ifelse(dl == 'Gaussian', ' with varying observation variances. ',
                          '. '),
                   'Colors represent different cluster (i.i.d nugget) variances',
                   ' used in an experiment.')

    if(mesh.resolution == 'all'){
      str1 <- paste(str1,
                    'Each point is the median bias of three experiments (low, medium, and',
                    'high resolution SPDE triangulations), calculated across 75 replicates,'
                    )
    }else{
      str1 <- paste0(str1,
                     ' Each point is the median bias of one experiment (with ', mesh.resolution,
                     ' resolution SPDE triangulation), calculated across 25 replicates,'
                     )
    }

    str1 <- paste(str1, 'and the bars represent the middle', cov.perc,
                  'quantile range of the bias across replicates.',
                  'Slight jitter has been applied to improve legibility.' )

    HTML(paste(title, str1, sep = '<br/>'))
  })

}

# make the app
shinyApp(ui = ui, server = server)
