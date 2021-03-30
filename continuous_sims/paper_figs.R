## plots made for TMB-INLA paper
## AOZ - 2020
require(raster)
require(viridis)

## ~~~~~~~~~~~~~~~~
## varying meshes ~
## ~~~~~~~~~~~~~~~~
source("./make_mesh_plot.R")

## ~~~~~~~~~~~~~~~~
## make cov plots ~
## ~~~~~~~~~~~~~~~~
cov_list <- readRDS('/media/FantomHD/IHME_DUMP/from_the_cluster/tmb_inla_sim/2020_08_08_13_13_30/cov_list_02.rds')

png(file = '~/Dropbox/AaronJon/TMB_SpatialStats/figures/cov_plot.png',
    width = 12, height = 5, units = 'in', res = 350)
par(mfrow = c(1, 2))
plot(cov_list[[1]], main = 'Access Time to Health Care',
     maxpixels = 1e9, col=rev(viridis(100)),
     legend.args=list(text='Hours', side=3, font=16, line=0, cex=.6))
plot(cov_list[[2]], main = '(Center-scaled) Malaria Incidence',
     maxpixels = 1e9, col=rev(magma(100)),
     legend.args=list(text='', side=3, font=16, line=0, cex = .8))
dev.off()
