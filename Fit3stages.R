##################################################
##                                              ##
## Fit the 3-stages model in FishExtremeCatcher ##
##           Using INLA_22.05.07                ##
##################################################
library(INLA)
library(rgdal)
load('FECdata.Rdata')
bathy = readGDAL("bathy.tiff")
bathy$band1[which(bathy$band1==0)] = NA

## Locs
spat.scale  = 1e05
locs        = matrix(c(FECdata$x/spat.scale,FECdata$y/spat.scale), ncol = 2)
## Mesh (provided by Iosu)
boundary = inla.nonconvex.hull(points = locs, convex = -.04)
load("mesh_extreme_fine.rdata")
mesh.points = SpatialPoints(coords=mesh$loc[,1:2]*spat.scale, proj4string=CRS(proj4string(bathy)))

##############################################
## Stage 1: intercept + bathy + year + spat ##
##############################################
## Spatial effect
spat.index = inla.spde.make.A(mesh, loc = locs)
spat.pred  = inla.spde.make.A(mesh)
spde_npue  = inla.spde2.pcmatern(mesh, prior.range=c(1.5, 0.15), prior.sigma=c(1, 0.2)) ### el mesh está en cientos de kilometros (spat.scale = 100 000), por lo tanto 1.5 es = 150km

## 1D SPDE for bathy
mesh1d    = inla.mesh.1d(seq(min(FECdata$profundidad), max(FECdata$profundidad), by = 10)) 
A1        = inla.spde.make.A(mesh1d, FECdata$profundidad)
spde1_NA  = inla.spde2.pcmatern(mesh1d, prior.range = c(150, NA), prior.sigma = c(1, 0.2))
spde1.idx = inla.spde.make.index("bathy", n.spde = spde1_NA$n.spde)

stack_stage1 = inla.stack(data = list(y=FECdata$numbers),
                          A = list(spat.index, 1,A1),
                          effects = list(spat = 1:mesh$n,
                                       list(intercept = 1,
                                            year = FECdata$yearID,
                                            effort = FECdata$inicio.virado),
                                       spde1.idx),
                          tag='est')
form1 = y ~ - 1 + intercept + f(bathy, model = spde1_NA)+ f(year, model = "rw2") + f(spat, model = spde_npue) + offset(log(effort))

stage1_fit = inla(form1,
                  family            = c("zeroinflatednbinomial2"),
                  data              = inla.stack.data(stack_stage1),
                  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
                  control.predictor = list(A = inla.stack.A(stack_stage1), compute = TRUE, link = 1),
                  verbose           = F, 
                  num.threads       = 1)

save(stage1_fit, mesh1d, file = "fitstage1.Rdata")
################################
## Stage 2: intercept + bathy ##
################################
fitted  = inla.stack.index(stack_stage1,"est")$data
means   = stage1_fit$summary.fitted.values[fitted,1]
perc_85 = unlist(lapply(means, function(x){qnbinom(p = 0.85,size = stage1_fit$summary.hyperpar[1,1], mu = x)}))

FECdata$extremes_presence85 = as.numeric(FECdata$numbers>=perc_85 & FECdata$npue > mean(FECdata$npue))

stack_stage2 = inla.stack(data = list(y = FECdata$extremes_presence85),
                          A    = list(spat.index, 1,A1),
                          effects = list(spat.bin = 1:mesh$n,
                                         list(intercept = 1,
                                              prof  = FECdata$profundidad), 
                                         spde1.idx),
                          tag='est')
form2 =  y ~ -1 + intercept +  prof 

stage2_fit = inla(form2,
                  family          = c('binomial'),
                  data            = inla.stack.data(stack_stage2), 
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                  control.predictor = list(A = inla.stack.A(stack_stage2), compute = TRUE),
                  verbose = F,
                  num.threads = 1)

save(stage2_fit, file = "fitstage2.Rdata")
########################################
## Stage 3: intercept +  year + bathy ##
########################################
exc = which(FECdata$extremes_presence85 == 1)
extreme_data = FECdata[exc, ]

## 1D SPDE for bathy
mesh1d_ext = inla.mesh.1d(seq(min(extreme_data$profundidad), max(extreme_data$profundidad), by = 10)) 
pred_bathy = over(mesh.points,bathy)$band1
pred_bathy[which(pred_bathy>max(FECdata$profundidad))] = max(FECdata$profundidad)
pred_bathy = pred_bathy
pred_bathy[which(is.na(pred_bathy))] = max(FECdata$profundidad)
pred_bathy[which(pred_bathy>max(extreme_data$profundidad))] = max(extreme_data$profundidad)

A1_ext        = inla.spde.make.A(mesh1d_ext, extreme_data$profundidad)
spde1_ext     = inla.spde2.pcmatern(mesh1d_ext, prior.range = c(80, 0.15), prior.sigma = c(2, 0.2))
spde1_NA_ext  = inla.spde2.pcmatern(mesh1d_ext, prior.range = c(80, NA), prior.sigma = c(2, 0.2))
spde1.ext.idx = inla.spde.make.index("bathy", n.spde = spde1_ext$n.spde)

i.exc = inla.spde.make.A(mesh, loc = locs[exc,])

stack_stage3 = inla.stack(data = list(y=extreme_data$numbers),
                          A = list(i.exc, 1,A1_ext),
                          effects = list(spat = 1:mesh$n,
                                       list(intercept = 1,
                                            year = extreme_data$yearID,
                                            effort = extreme_data$inicio.virado),
                                       spde1.ext.idx),
                          tag='est')

beta = 0.5 ## median
form3 =  y ~ -1 + intercept +  year + f(bathy, model = spde1_NA_ext) + offset(log(effort)) 

stage3_fit = inla(form3,
                  family            = c('dgp'),
                  control.family    = list(control.link = list(quantile=beta), hyper = list(tail = list(prior = "pc.gevtail", param = c(1, 0.0, 0.5)))),
                  data              = inla.stack.data(stack_stage3),
                  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
                  control.predictor = list(A=inla.stack.A(stack_stage3), compute=TRUE),
                  verbose           = F, 
                  num.threads       = 1)

save(stage3_fit, mesh1d_ext, extreme_data, file = "fitstage3.Rdata")
##########################################################################################
## Stage 3 with predictions at all observed locations (for predictive model assessment) ##
##########################################################################################
exc = which(FECdata$extremes_presence85 == 1)
extreme_data = FECdata
extreme_data$numbers[-exc] = NA

## 1D SPDE for bathy (same as stage1)
mesh1d    = inla.mesh.1d(seq(min(FECdata$profundidad), max(FECdata$profundidad), by = 10)) 
A1        = inla.spde.make.A(mesh1d, FECdata$profundidad)
spde1_NA  = inla.spde2.pcmatern(mesh1d, prior.range = c(150, NA), prior.sigma = c(1, 0.2))
spde1.idx = inla.spde.make.index("bathy", n.spde = spde1_NA$n.spde)

stack_stage3 = inla.stack(data = list(y=extreme_data$numbers),
                          A = list(spat.index, 1,A1),
                          effects = list(spat = 1:mesh$n,
                                         list(intercept = 1,
                                              prof = FECdata$profundidad,
                                              year = FECdata$yearID,
                                              effort = FECdata$inicio.virado),
                                         spde1.idx),
                          tag='est')

beta = 0.5 ## median
form3 =  y ~ -1 + intercept +  year + f(bathy, model = spde1_NA) + offset(log(effort)) 

stage3_pred = inla(form3,
                   family            = c('dgp'),
                   control.family    = list(control.link = list(quantile=beta), hyper = list(tail = list(prior = "pc.gevtail", param = c(1, 0.0, 0.5)))),
                   data              = inla.stack.data(stack_stage3),
                   control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
                   control.predictor = list(A=inla.stack.A(stack_stage3), compute=TRUE),
                   verbose           = F, 
                   num.threads       = 1)

save(stage3_pred, extreme_data, file = "predstage3.Rdata")



