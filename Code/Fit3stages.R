library(INLA)
library(rgdal)
load('Data/datos.Rdata')
bathy = readGDAL("Data/bathy.tiff")
bathy$band1[which(bathy$band1==0)] = NA

## Locs
spat.scale  = 1e05
locs        = matrix(c(datos$x/spat.scale,datos$y/spat.scale), ncol = 2)
## Mesh (provided by Iosu)
boundary = inla.nonconvex.hull(points = locs, convex = -.04)
load("Data/mesh_extreme_fine.rdata")
mesh.points = SpatialPoints(coords=mesh$loc[,1:2]*spat.scale, proj4string=CRS(proj4string(bathy)))

# ## Prediction 1D SPDE for bathymetry (stage 3)
mesh1d     = inla.mesh.1d(seq(min(datos$PROFUNDIDAD), max(datos$PROFUNDIDAD), by = 10))
pred_bathy = over(mesh.points,bathy)$band1
pred_bathy[which(pred_bathy>max(datos$PROFUNDIDAD))] = max(datos$PROFUNDIDAD)
pred_bathy_NA = pred_bathy
pred_bathy_NA[which(is.na(pred_bathy_NA))] = max(datos$PROFUNDIDAD)
# A1_pred = inla.spde.make.A(mesh1d,pred_bathy_NA)

## Prediction area (inside the polygon defined by locations)
bound    = Polygon(boundary$loc*spat.scale)
poly     = SpatialPolygons(list(Polygons(list(bound),'0')))
ea       = point.in.polygon(coordinates(mesh.points)[,1], coordinates(mesh.points)[,2], coordinates(bound)[,1], coordinates(bound)[,2], mode.checked = FALSE)
idx.mesh = which(ea%in%c(0))
plot(mesh$loc[-idx.mesh,1:2]*spat.scale)


##############################################
## Stage 1: intercept + bathy + year + spat ##
##############################################
## Spatial effect
spat.index = inla.spde.make.A(mesh, loc = locs)
spat.pred  = inla.spde.make.A(mesh)
spde_npue  = inla.spde2.pcmatern(mesh, prior.range=c(1.5, 0.15), prior.sigma=c(1, 0.2)) ### el mesh está en cientos de kilometros (spat.scale = 100 000), por lo tanto 1.5 es = 150km

## 1D SPDE for bathy
mesh1d    = inla.mesh.1d(seq(min(datos$PROFUNDIDAD), max(datos$PROFUNDIDAD), by = 10)) 
A1        = inla.spde.make.A(mesh1d, datos$PROFUNDIDAD)
spde1_NA  = inla.spde2.pcmatern(mesh1d, prior.range = c(150, NA), prior.sigma = c(1, 0.2))
spde1.idx = inla.spde.make.index("bathy", n.spde = spde1_NA$n.spde)

stack_stage1 = inla.stack(data = list(y=datos$numbers),
                          A = list(spat.index, 1,A1),
                          effects = list(spat = 1:mesh$n,
                                       list(intercept = 1,
                                            prof = datos$PROFUNDIDAD,
                                            year = datos$year2,
                                            effort = datos$INICIO.VIRADO),
                                       spde1.idx),
                          tag='est')
form1 = y ~ - 1 + intercept + f(bathy, model = spde1_NA)+ f(year, model = "rw2") + f(spat, model = spde_npue) + offset(log(effort))

stage1_fit = inla(form1,
                  family            = c("zeroinflatednbinomial2"),
                  data              = inla.stack.data(stack_stage1), 
                  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
                  control.predictor = list(A = inla.stack.A(stack_stage1), compute = TRUE, link = 1),
                  verbose           = F, 
                  num.threads       = 1, 
                  control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"))

save(stage1_fit, file = "Fits/fitstage1.Rdata")
################################
## Stage 2: intercept + bathy ##
################################
fitted  = inla.stack.index(stack_stage1,"est")$data
means   = stage1_fit$summary.fitted.values[fitted,1]
perc_85 = unlist(lapply(means, function(x){qnbinom(p = 0.85,size = stage1_fit$summary.hyperpar[1,1], mu = x)}))

datos$extremes_presence85 = as.numeric(datos$numbers>=perc_85 & datos$npue > mean(datos$npue))

stack_stage2 = inla.stack(data = list(y = datos$extremes_presence85),
                          A    = list(spat.index, 1,A1),
                          effects = list(spat.bin = 1:mesh$n,
                                       list(intercept = 1,
                                            prof  = datos$PROFUNDIDAD,
                                            prof2 = datos$PROFUNDIDAD^2,
                                            year  = datos$year2), ## change 
                                       spde1.idx),
                          tag='est')
form2 =  y ~ -1 + intercept +  prof 

stage2_fit = inla(form2,
                  family          = c('binomial'),
                  control.inla    = list(strategy = "simplified.laplace", int.strategy = "eb"),
                  data            = inla.stack.data(stack_stage2), 
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                  control.predictor = list(A = inla.stack.A(stack_stage2), compute = TRUE),
                  verbose = F,
                  num.threads = 1)
save(stage2_fit, file = "Fits/fitstage2.Rdata")
########################################
## Stage 3: intercept +  year + bathy ##
########################################
exc = which(datos$extremes_presence85 == 1)
extreme_data = datos[exc, ]
extreme_data$bathy_stage3 = inla.group(extreme_data$PROFUNDIDAD, n = 15, method = "cut")

## 1D SPDE for bathy
mesh1d_ext = inla.mesh.1d(seq(min(extreme_data$PROFUNDIDAD), max(extreme_data$PROFUNDIDAD), by = 10)) 

A1_ext        = inla.spde.make.A(mesh1d_ext, extreme_data$PROFUNDIDAD)
spde1_ext     = inla.spde2.pcmatern(mesh1d_ext, prior.range = c(80, 0.15), prior.sigma = c(2, 0.2))
spde1_NA_ext  = inla.spde2.pcmatern(mesh1d_ext, prior.range = c(80, NA), prior.sigma = c(2, 0.2))
spde1.ext.idx = inla.spde.make.index("bathy", n.spde = spde1_ext$n.spde)

pred_bathy_ext = pred_bathy_NA
pred_bathy_ext[which(pred_bathy_ext>max(extreme_data$PROFUNDIDAD))] = max(extreme_data$PROFUNDIDAD)
A1_pred_ext = inla.spde.make.A(mesh1d_ext,pred_bathy_ext)

i.exc = inla.spde.make.A(mesh, loc = locs[exc,])

stack_stage3 = inla.stack(data = list(y=extreme_data$numbers),
                          A = list(i.exc, 1,A1_ext),
                          effects = list(spat = 1:mesh$n,
                                       list(intercept = 1,
                                            year = extreme_data$year2,
                                            #month=datos$MonthShot[exc],
                                            prof = extreme_data$bathy_stage3,
                                            effort = extreme_data$INICIO.VIRADO),
                                       spde1.ext.idx),
                          tag='est')

alfa = 0.5 ## median
form3 =  y ~ -1 + intercept +  year + f(bathy, model = spde1_NA_ext) + offset(log(effort)) 

stage3_fit=inla(form3,
                family            = c('dgp'),
                control.inla      = list(strategy = "laplace"),
                # control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb"),
                control.family    = list(control.link = list(quantile=alfa), hyper = list(tail = list(prior = "pc.gevtail", param = c(1, 0.0, 0.5)))),
                data              = inla.stack.data(stack_stage3),
                control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
                control.predictor = list(A=inla.stack.A(stack_stage3), compute=TRUE),
                verbose           = F, 
                num.threads       = 1)
save(stage3_fit, extreme_data, file = "Fits/fitstage3.Rdata")
