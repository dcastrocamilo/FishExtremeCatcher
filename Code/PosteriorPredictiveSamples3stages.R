library(INLA)
library(rgdal)
library(ggplot2)
load('Data/datos.Rdata')
bathy = readGDAL("Data/bathy.tiff")
bathy$band1[which(bathy$band1==0)] = NA
load("Fits/fitstage1.Rdata")
load("Fits/fitstage2.Rdata")
load("Fits/fitstage3.Rdata")
load("Data/mesh_extreme_fine.rdata")

N = 500

## Prediction: 1D SPDE for bathymetry
mesh1d      = inla.mesh.1d(seq(min(datos$PROFUNDIDAD), max(datos$PROFUNDIDAD), by = 10)) 
spat.scale  = 1e05
locs        = matrix(c(datos$x/spat.scale,datos$y/spat.scale), ncol = 2)
mesh.points = SpatialPoints(coords=mesh$loc[,1:2]*spat.scale, proj4string=CRS(proj4string(bathy)))
pred_bathy = over(mesh.points,bathy)$band1
pred_bathy[which(pred_bathy>max(datos$PROFUNDIDAD))] = max(datos$PROFUNDIDAD)
pred_bathy_NA = pred_bathy
pred_bathy_NA[which(is.na(pred_bathy_NA))] = max(datos$PROFUNDIDAD)
# A1_pred = inla.spde.make.A(mesh1d,pred_bathy_NA)


##############################################
## Stage 1: intercept + bathy + year + spat ##
##############################################
psamp_stage1 = inla.posterior.sample(N, stage1_fit)
contents_stage1 = stage1_fit$misc$configs$contents

## Retain samples of the intercept
idx_Intercept = which(contents_stage1$tag == "intercept")
idx_Intercept = contents_stage1$start[idx_Intercept]:(contents_stage1$start[idx_Intercept] + contents_stage1$length[idx_Intercept]-1)
Intercept     = lapply(psamp_stage1, function(x)x$latent[idx_Intercept])
Intercept     = matrix(unlist(Intercept), ncol = 1, byrow = T)

## Retain samples of bathy
idx_bathy = which(contents_stage1$tag == "bathy")
idx_bathy = contents_stage1$start[idx_bathy]:(contents_stage1$start[idx_bathy]+contents_stage1$length[idx_bathy]-1)
bathy     = lapply(psamp_stage1, function(x) x$latent[idx_bathy])
bathy     = matrix(unlist(bathy), ncol = length(idx_bathy), byrow = T)

## Retain samples of year
idx_year = which(contents_stage1$tag == "year")
idx_year = contents_stage1$start[idx_year]:(contents_stage1$start[idx_year]+contents_stage1$length[idx_year]-1)
year     = lapply(psamp_stage1, function(x) x$latent[idx_year])
year     = matrix(unlist(year), ncol = length(idx_year), byrow = T)

## Retain samples of spatial effect
idx_i = which(contents_stage1$tag == "spat")
idx_i = contents_stage1$start[idx_i]:(contents_stage1$start[idx_i]+contents_stage1$length[idx_i]-1)
spat  = lapply(psamp_stage1, function(x) x$latent[idx_i])
spat  = matrix(unlist(spat), ncol = length(idx_i), byrow = T)

## Predicting at the mesh nodes
bathy_levels_stage1 = seq(min(datos$PROFUNDIDAD), max(datos$PROFUNDIDAD), by = 10)
year_levels_stage1 = stage1_fit$summary.random$year$ID

y = 1
mu1 = eta1 = matrix(NA, nrow = N, ncol = mesh$n)
prob0 = mean_mu1 = c()
for(i in 1:mesh$n){
  if(is.na(pred_bathy[i])){ # DATOS
    eta1[,i] = NA
    mu1[,i]  = NA
    prob0[i] = NA
  }else{
    idx.bathy = which.min(abs(bathy_levels_stage1-pred_bathy[i]))
    idx.year = y
    
    eta1[,i] = Intercept + bathy[,idx.bathy] + year[,idx.year] + spat[,i]
    mu1[,i]  = exp(eta1[,i])
    
    p0 = function(mu,alpha){1-(mu/(1+mu))^alpha}
    alpha=stage1_fit$summary.hyperpar[2,1]
    mean_mu1[i] = mean(mu1[,i])
    #mu = mean(exp(inla.rmarginal(n,numbers_pred$marginals.fitted.values[[pred.idx[i]]])))
    prob0[i] = 1 - p0(mu=mean_mu1[i],
                      alpha=alpha)
    
  }
}

p_data1 = data.frame(x = mesh$loc[,1], y = mesh$loc[,2], median_mu1 = apply(mu1,2,median),
                    mean_mu1 = apply(mu1,2,mean), prob = prob0)

ggplot(p_data1) + geom_point(aes(x=x,y=y,color=mean_mu1))
ggplot(p_data1) + geom_point(aes(x=x,y=y,color=median_mu1))
ggplot(p_data1) + geom_point(aes(x=x,y=y,color=prob))


################################
## Stage 2: intercept + bathy ##
################################
psamp_stage2    = inla.posterior.sample(N, stage2_fit)
contents_stage2 = stage2_fit$misc$configs$contents

## Retain samples of the intercept
idx_Intercept = which(contents_stage2$tag == "intercept")
idx_Intercept = contents_stage2$start[idx_Intercept]:(contents_stage2$start[idx_Intercept]+contents_stage2$length[idx_Intercept]-1)
Intercept     = lapply(psamp_stage2,function(x)x$latent[idx_Intercept])
Intercept     = matrix(unlist(Intercept),ncol=1,byrow=T)

## Retain samples of bathy
idx_bathy = which(contents_stage2$tag == "prof")
idx_bathy = contents_stage2$start[idx_bathy]:(contents_stage2$start[idx_bathy]+contents_stage2$length[idx_bathy]-1)
bathy     = lapply(psamp_stage2, function(x) x$latent[idx_bathy])
bathy     = matrix(unlist(bathy), ncol = length(idx_bathy), byrow = T)

## Predicting at the mesh nodes
inv.logit = function(x){exp(x)/(1+exp(x))}
y = 1
pi = prob = matrix(NA, nrow = N, ncol  =mesh$n)
mean_prob = c()
for(i in 1:mesh$n){
  if(is.na(pred_bathy[i])){
    pi[,i] = NA
    prob[,i] = NA
  }else{
    pi[,i] = Intercept + bathy*pred_bathy[i]
    prob[,i] = inv.logit(pi[,i])
    #mean_prob3[i] = mean(prob3[,i])
    #mu = mean(exp(inla.rmarginal(n,numbers_pred$marginals.fitted.values[[pred.idx[i]]])))
  }
}

p_data2 = data.frame(x = mesh$loc[,1],y = mesh$loc[,2], mean_prob = apply(prob,2,mean))

ggplot(p_data2) + geom_point(aes(x=x,y=y,color=mean_prob))


########################################
## Stage 3: intercept +  year + bathy ##
########################################
psamp_stage3    = inla.posterior.sample(N, stage3_fit)
contents_stage3 = stage3_fit$misc$configs$contents

## Retain samples of the intercept
idx_Intercept = which(contents_stage3$tag == "intercept")
idx_Intercept = contents_stage3$start[idx_Intercept]:(contents_stage3$start[idx_Intercept]+contents_stage3$length[idx_Intercept]-1)
Intercept     = lapply(psamp_stage3,function(x)x$latent[idx_Intercept])
Intercept     = matrix(unlist(Intercept),ncol=1,byrow=T)

## Retain samples of bathy
idx_bathy = which(contents_stage3$tag == "bathy")
idx_bathy = contents_stage3$start[idx_bathy]:(contents_stage3$start[idx_bathy]+contents_stage3$length[idx_bathy]-1)
bathy     = lapply(psamp_stage3, function(x) x$latent[idx_bathy])
bathy     = matrix(unlist(bathy), ncol = length(idx_bathy), byrow = T)

## Retain samples of year
idx_year = which(contents_stage3$tag == "year")
idx_year = contents_stage3$start[idx_year]:(contents_stage3$start[idx_year]+contents_stage3$length[idx_year]-1)
year     = lapply(psamp_stage3, function(x) x$latent[idx_year])
year     = matrix(unlist(year), ncol = length(idx_year), byrow = T)

## Predicting at the mesh nodes
bathy_levels_stage3 = seq(min(extreme_data$PROFUNDIDAD), max(extreme_data$PROFUNDIDAD), by = 10)
i = 1
y = 1
mu2 = eta2 = matrix(NA,nrow = N, ncol = mesh$n)
mean_mu2 = c()
for(i in 1:mesh$n){
  if(is.na(pred_bathy[i])){
    eta2[,i] = NA
    mu2[,i] = NA
  }else{
    idx.bathy = which.min(abs(bathy_levels_stage3-pred_bathy[i]))
    idx.year = y
    
    eta2[,i] = Intercept + bathy[,idx.bathy] + year*idx.year
    mu2[,i] = exp(eta2[,i])
    
    #mean_mu2[i] = mean(mu2[,i])
    
  }
}

p_data3 = data.frame(x = mesh$loc[,1],y = mesh$loc[,2], median_mu2 = apply(mu2,2,median),
                     mean_mu2 = apply(mu2, 2, mean))

ggplot(p_data3) + geom_point(aes(x=x,y=y,color=mean_mu2))


