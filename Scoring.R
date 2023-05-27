##########################################################################
##                                                                      ##
## Predictive model assessment for the baseline and the 3 stages models ##
##                                                                      ##
##########################################################################
library(gamlss.dist) # for rZINBI. Note that sigma = 1/n_NB and nu = 1 - (mu/(1+mu))^k_NB
library(INLA)
library(kdensity)
library(rje)
load('FECdata.Rdata')
load("fitstage1.Rdata")
load("fitstage2.Rdata")
load("predstage3.Rdata")
y = FECdata$numbers
N = 500
beta = 0.5

#########################
## Auxiliary functions ##
#########################
# obtain each parameter/hyperparameter from inla.posterior.sample()
get.posterior.param.baseline = function(psamp, idx2keep, which = 'all'){
  mu_NB = exp(psamp$latent[idx2keep, 1])
  n_NB  = psamp$hyperpar['size for nbinomial zero-inflated observations']
  k_NB  = psamp$hyperpar['parameter alpha for zero-inflated nbinomial2']
  pr0   = 1 - (mu_NB/(1+mu_NB))^k_NB
  
  mean_ZINB = (1-pr0)*mu_NB
  var_ZINB  = (1-pr0)*mean_ZINB*(1 + mean_ZINB*(1/n_NB + pr0))
  
  if(which == 'mean_ZINB')
    return(as.numeric(mean_ZINB))
  if(which == 'var_ZINB')
    return(as.numeric(var_ZINB))
  if(which == 'sd_ZINB')
    return(as.numeric(sqrt(var_ZINB)))
  if(which == 'mu')
    return(as.numeric(mu_NB))
  if(which == 'sigma')
    return(as.numeric(1/n_NB))
  if(which == 'pr0')
    return(as.numeric(pr0))
  if(which == 'all')
    return(list(mean_ZINB = as.numeric(mean_ZINB), var_ZINB = as.numeric(var_ZINB), sd_ZINB = as.numeric(sqrt(var_ZINB)),
                mu_NB = as.numeric(mu_NB), sigma_NB = as.numeric(1/n_NB), pr0 = as.numeric(pr0)))
}

# samples from baseline predictive posterior distribution (for data)
sim.from.ppdistribution.baseline = function(psamp, idx2keep) {
  params = get.posterior.param.baseline(psamp, idx2keep)
  
  y = rep(NA, length(idx2keep))
  for(i in 1:length(idx2keep)){
    # in rZINBI, mu and sigma are parameters of the NBI
    y[i] = rZINBI(1, mu = params$mu_NB[i], sigma = params$sigma_NB, nu = params$pr0[i])
  }
  y
}

# random generation of DGPD data
rdpgd = function(n, xi, l){
  if(any(xi <= 0, l <= 0)){
    stop('xi must be positive')
  }else{
    u = runif(n)
    y = (1/(1-u)^xi-1)*l/xi-1
    as.integer(y)
  }
}

# samples from 3stages predictive posterior distribution (for data)
sim.from.ppdistribution.3stages = function(psamp1, psamp2, psamp3, idx2keep) {
  
  # Params stage1
  mu_NB = exp(psamp1$latent[idx2keep, 1])
  n_NB  = psamp1$hyperpar['size for nbinomial zero-inflated observations']
  k_NB  = psamp1$hyperpar['parameter alpha for zero-inflated nbinomial2']
  pr0   = 1 - (mu_NB/(1+mu_NB))^k_NB
  
  # Params stage2
  q_Be = exp(psamp2$latent[idx2keep, 1])/(1+exp(psamp2$latent[idx2keep, 1]))
  
  # Params stage3
  mBeta_DGP = exp(psamp3$latent[idx2keep, 1])
  xi_DGP    = psamp3$hyperpar['Tail parameter for the dgp observations']
  l_DGP     = xi_DGP*mBeta_DGP/(1/(1-beta)^xi_DGP - 1)
  
  
  y = rep(NA, length(mu_NB))
  for(i in 1:length(mu_NB)){
    z = rbinom(1, size = 1, prob = q_Be[i])
    if(z == 0){
      y[i] = rZINBI(1, mu = mu_NB[i], sigma = 1/n_NB, nu = pr0[i])  
    }
    if(z == 1){
      y[i] = rdpgd(1, xi = xi_DGP, l = l_DGP[i])
    }
  }
  y
}

#############################################################
## Posteriors from the linear predictors & hyperparameters ##
#############################################################
psamp_stage1    = inla.posterior.sample(N, stage1_fit)
psamp_stage2    = inla.posterior.sample(N, stage2_fit)
psamp_stage3    = inla.posterior.sample(N, stage3_pred)

#################################
## Predictive model assessment ##
#################################
contents_stage1 = stage1_fit$misc$configs$contents
idx2keep        = 1:contents_stage1$length[1]
Psamp.baseline  = sapply(psamp_stage1, sim.from.ppdistribution.baseline, idx2keep = idx2keep)
Psamp.3stages   = matrix(NA, length(y), N)
for(j in 1:N){
  Psamp.3stages[,j] = sim.from.ppdistribution.3stages(psamp_stage1[[j]], psamp_stage2[[j]], psamp_stage3[[j]], idx2keep = idx2keep)
}

scoringBaseline = function(y, Psamp){
  # See https://cran.r-project.org/web/packages/tscount/tscount.pdf and 
  # Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for count data. Biometrics 65, 1254–1261
  
  #- y (vector, length n) observations
  #- Psamp (matrix, n x N) i-th row contains a random sample of size N of the posterior predictive distribution of y[i]
  
  Pdist = Pdens = rep(NA, length(y))
  x = 0:5e03
  Px = px = ind = matrix(NA, length(y), length(x))
  for(i in 1:length(y)){ # takes a while...
    printPercentage(i, length(y))
    eCDF     = ecdf(Psamp.baseline[i,])
    kden     = kdensity(Psamp.baseline[i,])
    Pdist[i] = eCDF(y[i])
    Pdens[i] = kden(y[i])
    Px[i,]   = eCDF(x)
    px[i,]   = kden(x)
    ind[i,]  = as.numeric(y[i] <= x)
  }
  p2   = sum(px^2)
  
  logs = -mean( log(Pdens.baseline) ) # average log-score
  qs   = mean( -2*Pdens.baseline + p2 ) # average quadratic score
  sphs = mean( -Pdens.baseline/sqrt(p2) ) # average spherical score
  rps  = mean( rowSums((Px - ind)^2) )
  list(logs = logs, qs = qs, sphs = sphs, rps = rps)
}

scoring3stages = function(y, Psamp){
  #- y (vector, length n) observations
  #- Psamp (matrix, n x N) i-th row contains a random sample of size N of the posterior predictive distribution of y[i]
  
  Pdist = Pdens = rep(NA, length(y))
  x = 0:5e03
  Px = px = ind = matrix(NA, length(y), length(x))
  for(i in 1:length(y)){ # takes a while...
    printPercentage(i, length(y))
    eCDF     = ecdf(Psamp[i,])
    kden     = tryCatch(kdensity(Psamp[i,]), error = function(e) e)
    if(inherits(kden, "error"))
      kden     = kdensity(Psamp[i,], kernel = "gamma")
    Pdist[i] = eCDF(y[i])
    Pdens[i] = kden(y[i])
    Px[i,]   = eCDF(x)
    px[i,]   = kden(x)
    ind[i,]  = as.numeric(y[i] <= x)
  }
  p2   = sum(px^2)
  
  logs = -mean( log(Pdens) ) # average log-score
  qs   = mean( -2*Pdens + p2 ) # average quadratic score
  sphs = mean( -Pdens/sqrt(p2) ) # average spherical score
  rps  = mean( rowSums((Px - ind)^2) )
  list(logs = logs, qs = qs, sphs = sphs, rps = rps)
  
}

# Run
system.time(Sbaseline <- scoringBaseline(y, Psamp.baseline))
system.time(S3stages <- scoring3stages(y, Psamp.3stages))
Sbaseline
S3stages



