## this file preps objects meant for continuout sims for the discrete sims
## aoz 2020oct23

setwd("~/Documents/GitRepos/tmb_inla_comp")

nga.data.dir <- '/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects/'

# load pkgs

require('raster')
require('sp')
require('glue')
require('data.table')
require('spdep')


# load gaul_list, pop_raster, simple_polygon, simple_raster, subset_shape
load(file.path(nga.data.dir, 'poly_shape.rdata'))
pop_raster <- pop_raster[[1]]

# load nga admin1
nga1 <- shapefile(file.path(nga.data.dir, 'nga_shp/nga_admbnda_adm1_osgof_20190417.shp'))

# load covs
cov.list <- readRDS('/home/merms/Documents/Research/2020_tmb_v_inla/nigeria_objects/cov_list_02.rds')

# rasterize nga1, after making integer codes
nga1$ADM1_CODE <- as.numeric(gsub(pattern = 'NG', replacement = '', nga1$ADM1_PCODE))
nga1.ras <- rasterize(nga1, simple_raster, field='ADM1_CODE')
# plot(nga1.ras)

# aggregated pop and make pop-wgted covs for each ad1
agg.cov.list <- rbindlist(lapply(nga1$ADM1_CODE, function(code){
  sum.pop <- sum(pop_raster[nga1.ras==code], na.rm=T)
  mean.access <- sum(cov.list[['access2']][nga1.ras==code] *
                       pop_raster[[1]][nga1.ras==code], na.rm=T) / sum.pop
  mean.mapincidence <- sum(cov.list[['mapincidence']][nga1.ras==code] *
                             pop_raster[nga1.ras==code], na.rm=T) / sum.pop
  data.frame('adm1_code'=code,
             'agg_pop'=sum.pop,
             'agg_access'=mean.access,
             'agg_mapinc'=mean.mapincidence)}))

# merge on to nga1 shp
nga1 <- merge(nga1, agg.cov.list, by.x='ADM1_CODE', by.y='adm1_code')
## head(nga1@data)
## brks <- quantile(nga1$agg_pop, seq(0,1,1/7))
## cols_gr <- grey((length(brks):2)/length(brks))
## cols <- viridis((length(brks) - 1))
## dens <- rev((2:length(brks))*5)
## plot(nga1, col=cols[findInterval(nga1$agg_pop, brks, all.inside=TRUE)])
## plot(nga1, density=dens[findInterval(nga1$agg_pop, brks, all.inside=TRUE)])
## plot(nga1, col = "grey", boundary = "black")

# get the adjacency matrix for the admin1 shapes
nga1.adj <- poly2nb(nga1)

# and get the binary and row-standardized adj mats
W.nga1    <- nb2mat(nga1.adj, style='B')
W.nga1.rs <- nb2mat(nga1.adj, style='W')

# save prepped objects
save(nga1, nga1.ras, nga1.adj, W.nga1, W.nga1.rs, file = file.path(nga.data.dir, "prepped_discrete_obj.RData"))
