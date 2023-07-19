#############################
### Compute DIC and WAIC
#############################

######################
### Model: sp
######################

rm(list=ls())
load("rainfall_sp.RData")

dim(data_train)
str(m.1)
head(m.1$p.theta.samples)
head(m.1$p.beta.samples)

library(fields)

beta_post <- m.1$p.beta.samples[1001:2000,]
sigma_post <- m.1$p.theta.samples[1001:2000,1]
tau_post <- m.1$p.theta.samples[1001:2000,2]
phi_post <- m.1$p.theta.samples[1001:2000,3]

library(mvtnorm)

######## DIC

llikm <- c()
#iterate for j (MCMC iterations)
for(j in 1:nrow(beta_post))
{
  #iterate for i (observation)
  llikn <- c()
  for(i in 1:nrow(data_train))
  {
    #compute log(p(y|theta))
    llikn[i] <- dnorm(data_train$data[i],mean=sum(beta_post[j,] * X_train[i,]),
                          sd=sigma_post[j],log=T)
  }
  llikn <- llikn[llikn != -Inf]
  llikm[j] <- sum(llikn)
}

#compute Dbar
Dbar_sp <- -2 * sum(llikm)/nrow(beta_post)

#compute D1 (deviance in the first definition of Celeux et al (2006))
mu_mean <- X_train %*% matrix(apply(beta_post,2,mean),ncol=1)
sig_mean <- mean(sigma_post)
#iterate for i (observation)
llikn <- c()
for(i in 1:nrow(data_train))
{
  #compute log(p(y|theta))
  llikn[i] <- dnorm(data_train$data[i],mean=mu_mean[i],
                        sd=sig_mean,log=T)
}
llikn <- llikn[llikn != -Inf]

D1_sp <- -2 * sum(llikn)

DIC_sp <- 2 * Dbar_sp - D1_sp

save.image("DIC_WAIC.RData")

######## WAIC

#compute LPPD

#iterate for i (observations)
LPPDvec_sp <- c()
for(i in 1:nrow(data_train))
{
  #iterate for j (MCMC iterations)
  llikm <- c()
  for(j in 1:nrow(beta_post))
  {
    #compute p(y|theta)
    llikm[j] <- dnorm(data_train$data[i],mean=sum(beta_post[j,] * X_train[i,]),
                              sd=sigma_post[j])
  }
  LPPDvec_sp[i] <- log(mean(llikm))
}
LPPDvec_sp <- LPPDvec_sp[LPPDvec_sp != -Inf]
LPPD_sp <- sum(LPPDvec_sp)
LPPD_sp

save.image("DIC_WAIC.RData")

#compute pWAIC

llikn <- c()
for(i in 1:nrow(data_train))
{
  likm <- c()
  llikm <- c()
  #iterate for m (MCMC samples)
  for(j in 1:nrow(beta_post))
  {
    llikm[j] <- dnorm(data_train$data[i],mean=sum(beta_post[j,] * X_train[i,]),
                      sd=sig_mean,log=T)
    likm[j] <- exp(llikm[j]) 
  }
  llikn[i] <- log(mean(likm)) - mean(llikm[llikm != -Inf])
}  
llikn <- llikn[llikn != -Inf]
head(llikn)
pWAIC_sp <- 2 * sum(llikn[is.na(llikn) == F])

WAIC_sp <- 2 * pWAIC_sp - 2 * LPPD_sp

save.image("DIC_WAIC.RData")

######################
### Model: sSB
######################

load("rainfall_SB.RData")
str(obj_sSB)

mu_post <- obj_sSB$mu_post[(burn+1):runs,]
head(mu_post)

sig_post <- obj_sSB$sige_post[(burn+1):runs,]
head(sig_post)

pi_post <- obj_sSB$probs_last

lmix_dens <- function(y,mus,sigs,weights)
{
  mixd <- 0
  for(k in 1:length(mus)){
    mixd <- mixd + weights[k] * dnorm(y,mus[k],sigs[k])
  }
  return(log(mixd))
}

######## DIC

i <- 1
j <- 1

llikm <- c()
#iterate for j (MCMC iterations)
for(j in 1:nrow(mu_post))
{
  #iterate for i (observation)
  llikn <- c()
  for(i in 1:nrow(data_train))
  {
    #compute log(p(y|theta))
    llikn[i] <- lmix_dens(data_train$data[i],mus=mu_post[j,],
                          sigs=sig_post[j,],weights=pi_post[i,])
  }
  llikn <- llikn[llikn != -Inf]
  llikm[j] <- sum(llikn)
}

#compute Dbar
Dbar_sSB <- -2 * sum(llikm)/nrow(mu_post)

#compute D1 (deviance in the first definition of Celeux et al (2006))
mu_mean <- apply(mu_post,2,mean)
sig_mean <- apply(sig_post,2,mean)
#iterate for i (observation)
llikn <- c()
for(i in 1:nrow(data_train))
{
  #compute log(p(y|theta))
  llikn[i] <- lmix_dens(data_train$data[i],mus=mu_mean,
                        sigs=sig_mean,weights=pi_post[i,])
}
llikn <- llikn[llikn != -Inf]

D1_sSB <- -2 * sum(llikn)

DIC_sSB <- 2 * Dbar_sSB - D1_sSB

save.image("DIC_WAIC.RData")

######## WAIC

#compute LPPD

#iterate for i (observations)
LPPDvec_sSB <- c()
for(i in 1:nrow(data_train))
{
  #iterate for j (MCMC iterations)
  llikm <- c()
  for(j in 1:nrow(mu_post))
  {
    #compute p(y|theta)
    llikm[j] <- exp(lmix_dens(data_train$data[i],mus=mu_post[j,],
                              sigs=sig_post[j,],weights=pi_post[i,]))
  }
  LPPDvec_sSB[i] <- log(mean(llikm))
}
LPPDvec_sSB <- LPPDvec_sSB[LPPDvec_sSB != -Inf]
LPPD_sSB <- sum(LPPDvec_sSB)
LPPD_sSB

save.image("DIC_WAIC.RData")

#compute pWAIC

llikn <- c()
for(i in 1:nrow(data_train))
{
  likm <- c()
  llikm <- c()
  #iterate for m (MCMC samples)
  for(j in 1:nrow(mu_post))
  {
    llikm[j] <- lmix_dens(data_train$data[i],mus=mu_post[j,],
                              sigs=sig_post[j,],weights=pi_post[i,])
    likm[j] <- exp(llikm[j]) 
  }
  llikn[i] <- log(mean(likm)) - mean(llikm[llikm != -Inf])
}  
head(llikn)
pWAIC_sSB <- 2 * sum(llikn[is.na(llikn) == F])

WAIC_sSB <- 2 * pWAIC_sSB - 2 * LPPD_sSB

save.image("DIC_WAIC.RData")

######################
### Model: stSB
######################

load("rainfall_stSB.RData")
str(obj_stSB)

mu_post <- obj_stSB$mu_post[(burn+1):runs,]
head(mu_post)

sig_post <- obj_stSB$sige_post[(burn+1):runs,]
head(sig_post)

pi_post <- obj_stSB$probs_last

lmix_dens <- function(y,mus,sigs,weights)
{
  mixd <- 0
  for(k in 1:length(mus)){
    mixd <- mixd + weights[k] * dnorm(y,mus[k],sigs[k])
  }
  return(log(mixd))
}

######## DIC

llikm <- c()
#iterate for j (MCMC iterations)
for(j in 1:nrow(mu_post))
{
  #iterate for i (observation)
  llikn <- c()
  for(i in 1:nrow(data_train))
  {
    #compute log(p(y|theta))
    llikn[i] <- lmix_dens(data_train$data[i],mus=mu_post[j,],
                          sigs=sig_post[j,],weights=pi_post[i,])
  }
  llikn <- llikn[llikn != -Inf]
  llikm[j] <- sum(llikn)
}

#compute Dbar
Dbar_stSB <- -2 * sum(llikm)/nrow(mu_post)

#compute D1 (deviance in the first definition of Celeux et al (2006))
mu_mean <- apply(mu_post,2,mean)
sig_mean <- apply(sig_post,2,mean)
#iterate for i (observation)
llikn <- c()
for(i in 1:nrow(data_train))
{
  #compute log(p(y|theta))
  llikn[i] <- lmix_dens(data_train$data[i],mus=mu_mean,
                        sigs=sig_mean,weights=pi_post[i,])
}
llikn <- llikn[llikn != -Inf]

D1_stSB <- -2 * sum(llikn)

DIC_stSB <- 2 * Dbar_stSB - D1_stSB

save.image("DIC_WAIC.RData")

######## WAIC

#compute LPPD

#iterate for i (observations)
LPPDvec_stSB <- c()
for(i in 1:nrow(data_train))
{
  #iterate for j (MCMC iterations)
  llikm <- c()
  for(j in 1:nrow(mu_post))
  {
    #compute p(y|theta)
    llikm[j] <- exp(lmix_dens(data_train$data[i],mus=mu_post[j,],
                              sigs=sig_post[j,],weights=pi_post[i,]))
  }
  LPPDvec_stSB[i] <- log(mean(llikm))
}
LPPDvec_stSB <- LPPDvec_stSB[LPPDvec_stSB != -Inf]
LPPD_stSB <- sum(LPPDvec_stSB)
LPPD_stSB

save.image("DIC_WAIC.RData")

#compute pWAIC

llikn <- c()
for(i in 1:nrow(data_train))
{
  likm <- c()
  llikm <- c()
  #iterate for m (MCMC samples)
  for(j in 1:nrow(mu_post))
  {
    llikm[j] <- lmix_dens(data_train$data[i],mus=mu_post[j,],
                          sigs=sig_post[j,],weights=pi_post[i,])
    likm[j] <- exp(llikm[j]) 
  }
  llikn[i] <- log(mean(likm)) - mean(llikm[llikm != -Inf])
}  
head(llikn)
pWAIC_stSB <- 2 * sum(llikn[is.na(llikn) == F])

WAIC_stSB <- 2 * pWAIC_stSB - 2 * LPPD_stSB

save.image("DIC_WAIC.RData")

