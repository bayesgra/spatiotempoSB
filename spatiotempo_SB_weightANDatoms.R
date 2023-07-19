#############################################################################
# Do MCMC sampling for the spatial stick-breaking model
#    Here's the model:
#    
#    y[s] ~ N(x*beta+mu[g[s]],sige[g[s])   
#    mu[k]~dnorm(sigs), k=1,...n.terms
#    sige[1],...,sige[n.terms],sigs~dunif(0,mx.sig)
#
#    g[s] ~ Categorical(p_1[s],...,p_m[s])
#    p_j[s] = v_j[s] * prod_{k<j} (1-v_k[s])
#    v_k[s] = w_j(knots[j,],x[s,],rho) * v[k]
#            (w is a Gaussian kernel)
#
#    knot[j1,j2] ~ unif(0,1)
#    rho ~ Unif(0,1)
#    v[j] ~ beta(1,DOF)
#
#############################################################################

library(MASS)
library(SpatialTools)
library(mvtnorm)

#Betagn spike and slab prior
prior_beta <- function(betagn,p=rep(1/3,3),alpha=1,beta=1){
  p[1] * (betagn==0) + p[2] * (betagn==1) + p[3] * dbeta(betagn,alpha,beta)
}

#Gneting kernel
gneiting <- function(a,t,d,beta)
{
  1/(a*abs(t)+1) * exp(-d/((a*abs(t)+1)^(beta/2)))
}

#Compute v_k(s)
onevs<-function(z,rho,knot,v,kernel="gaussian",t=NULL,a=0,beta=0){
  kkk<-v
  for(j in 1:length(v)){
    if(kernel=="gaussian"){
      kkk<-kkk*rho^((z[,j]-knot[j])^2)
    }     
    if(kernel=="uniform"){
      kkk<-kkk*ifelse(abs(z[,j]-knot[j])<rho,1,0.01)
    }     
    if(kernel=="gneiting"){
      kkk <- kkk * gneiting(a,t,(z[,j]-knot[j])^2,beta)
    }
  }
  kkk}   

# n=50
# beta<-rep(0,p)
# v<-rep(.9,n.terms)
# sige<-rep(mx.sige/2,n.terms);taue<-1/sige^2
# mu<-rep(0,n.terms)
# sigs<-mx.sigs/2;taumu<-1/sigs^2
# knot<-matrix(runif(2*n.terms,0,1),2,n.terms)
# rho<-.5
# g<-rep(1,n)
# vs<-matrix(0,n,n.terms)
# 
# s1 <- seq(0,1,length=50)
# s2 <- seq(0,1,length=50)
# z <- cbind(sample(s1,50,replace=T), sample(s2,50,replace=T))
# t <- rep(1:10,5)
# 
# for(k in 1:n.terms){vs[,k]<-onevs(z,rho,knot[,k],v[k],kernel="gneiting",t=t,a=1,beta=1)}
# 
# probs<-makeprobs(vs)

#take all the vs and compute probabilities
makeprobs<-function(vs){
  m<-ncol(vs)
  failures<-matrix(0,m,m)
  for(j in 2:m){failures[1:(j-1),j]<-1}
  probs<-exp(log(vs)+log(1-vs)%*%failures)
  probs[,m]<-1-apply(probs[,1:m],1,sum)
  probs}

#generate categorical variables:
rcat<-function(prob){
  (1:length(prob))%*%rmultinom(1,1,ifelse(prob>0,prob,0))
}


#The whole Shebang!
STSB_varyingatoms<-function(y,x=NA,z,t=NULL,DOF=1,mx.sige=1,mx.sigs=1,n.terms=50,pred=F,x_new,z_new,t_new=NULL,
              runs=5000,burn=1000,display=10,kernel="gaussian",a=1,betagn=NULL){
  #y:       data
  #x:       covariates
  #z:       n x 2 matrix of coordinates, scaled to [0,1]
  #t:       vector of times
  #n.terms: number of terms in the mixture dist.
  #runs:    number of MCMC samples
  #burn:    number of samples to discard
  #display: how often to display results
  #DOF:     v~beta(1,DOF)
  #mx.sig:  the sds are U(0,mx.sd)
  #pred:    if there is prediction
  #x_new:   in case pred=T, covariates for prediction
  #z_new:   in case pred=T, coordinates for prediction
  #a:       parameter of the gneiting kernel, default at 1
  #betagn:  parameter of the gneiting kernel. If not given, it is meant to be an unknown parameter
  
  n<-length(y)
  if(is.na(max(x))){x<-matrix(1,n,1)}
  p<-ncol(x)
  
  #Standardize the outcomes and predictors:
  if(min(z)<0 | max(z)>1){print("THE SPATIAL COORDINATES ARE NOT STANDARDIZED!")}   
  
  #initial values
  beta<-rep(0,p)
  v<-rep(.9,n.terms)
  sige<-rep(mx.sige/2,n.terms);taue<-1/sige^2
  mu<-rep(0,n.terms)
  sigs<-mx.sigs/2;taumu<-1/sigs^2
  knot<-matrix(runif(2*n.terms,0,1),2,n.terms)
  rho<-.5
  g<-rep(1,n)
  vs<-matrix(0,n,n.terms)
  if(is.null(betagn)==T){
    betagn <- 0
  }
#  for(k in 1:n.terms){vs[,k]<-onevs(z,rho,knot[,k],v[k])}
  for(k in 1:n.terms){vs[,k]<-onevs(z,rho,knot[,k],v[k],kernel=kernel,t=t,a=1,
                                    beta=betagn)} 
  probs<-makeprobs(vs)
  
  sumtp<-summu<-summu2<-rep(0,n)
  count<-afterburn<-0
  keeprho<-rep(0,runs)
  keepbeta<-matrix(0,runs,p)
  keepbetagn <- rep(0,runs)
  
  mu_post <- matrix(0,runs,n)
  sige_post <- matrix(0,runs,n.terms)
  
  y_new <- matrix(0,runs,length(x_new))
  for(i in 1:runs){
    if(i/100 == round(i/100))(print(i))
    #Update beta
    COV<-solve(t(x)%*%diag(taue[g])%*%x)
    mn<-COV%*%t(x)%*%diag(taue[g])%*%(y-mu[g])
    beta<-mn+t(chol(COV))
    r<-y-x%*%beta
    
    mu <- rep(NA, n)
    for(j in 1:n.terms){
      #update mu
      nobs<-sum(g==j)
      mu[j]<-rnorm(1,0,sigs)
      
      if(nobs>1){
        varj <- 1/(taumu+nobs*taue[j])
        covj <- cov.st(coords=coords_train[g==j,], time=t[g==j], sp.type = "exponential", 
                       sp.par = c(varj,1))
        temp_mu <-mvrnorm(1,r[g==j],covj$V)
        mu[g==j] <- temp_mu}
      if(nobs==1){
        temp_mu <- rnorm(1,r[g==j],1/sqrt(taumu+nobs*taue[j]))
        mu[g==j] <- temp_mu
      }

      #update sige
      cansige<-sige;cansige[j]<-rnorm(1,sige[j],.1)
      if(sum(g==j)>1 & cansige[j]>0 & cansige[j]<mx.sige){
        cancovj <- cov.st(coords=coords_train[g==j,], time=t[g==j], sp.type = "exponential", 
                          sp.par = c(cansige[j],1))
        MHrate<-sum(dmvnorm(r[g==j],mu[g==j],cancovj$V,log=T)-
                      dmvnorm(r[g==j],mu[g==j],covj$V,log=T))
#        MHrate<-sum(dnorm(r[g==j],mu[j],cansige[j],log=T)-
#                      dnorm(r[g==j],mu[j],sige[j],log=T))
       if(runif(1,0,1)<exp(MHrate)){sige<-cansige} 
      }
      if(sum(g==j)<2 & cansige[j]>0 & cansige[j]<mx.sige){
        sige<-cansige 
      }
      taue<-1/sige^2
    }
    
    if(i>burn){
      mu_post[i,] <- mu
      sige_post[i,] <- sige
    }
    
    #update sigs:
    cansigs<-rnorm(1,sigs,.1)
    if(cansigs>0 & cansigs<mx.sigs){
      MHrate<-sum(dnorm(mu,0,cansigs,log=T)-
                    dnorm(mu,0,sigs,log=T))
      if(runif(1,0,1)<exp(MHrate)){sigs<-cansigs} 
    }
    taumu<-1/sigs^2
    
    #update g
    for(s in 1:n){
      cang<-g;cang[s]<-rcat(probs[s,])
      MHrate<-dnorm(r[s],mu[cang[s]],sige[cang[s]],log=T)-
        dnorm(r[s],mu[g[s]],sige[g[s]],log=T)
      if(runif(1,0,1)<exp(MHrate)){g<-cang} 
    }
    
    if(kernel=="gneiting"){
      #update betagn:
      ubeta <- runif(1,0,1)
      if(ubeta < 1/3){
        canbetagn <- 0
      } else {
        if(ubeta > 2/3){
          canbetagn <- runif(1,0,1)
        } else {
          canbetagn <- 1
        }
      }
      canvs<-vs 
      for(k in 1:n.terms){canvs[,k]<-onevs(z,rho,knot[,k],v[k],kernel="gneiting",t=t,a=1,
                                           beta=canbetagn)}
      canprobs<-makeprobs(canvs)
      MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                    log(probs[cbind(1:n,g)])) + log(prior_beta(canbetagn) - log(prior_beta(betagn)))
      if(runif(1,0,1)<exp(MHrate)){
        betagn<-canbetagn;vs<-canvs;probs<-canprobs}
      
      } else {
        #update rho:
        canrho<-runif(1,0,1);canvs<-vs 
        for(k in 1:n.terms){canvs[,k]<-onevs(z,canrho,knot[,k],v[k])}
        canprobs<-makeprobs(canvs)
        MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                      log(probs[cbind(1:n,g)]))     
        if(runif(1,0,1)<exp(MHrate)){
          rho<-canrho;vs<-canvs;probs<-canprobs}
      }
    
    #update v:
    for(k in 1:(n.terms-1)){
      if(max(g)<k){v[k]~rbeta(1,1,DOF)}
      if(max(g)>=k){
        canv<-v;canv[k]<-rnorm(1,v[k],.05)
        if(canv[k]>0 & canv[k]<1){
          canvs[,k]<-onevs(z,rho,knot[,k],canv[k],kernel=kernel,t=t,a=1,beta=betagn)
          canprobs<-makeprobs(canvs)
          MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                        log(probs[cbind(1:n,g)]))+
            log(dbeta(canv[k],1,DOF)/dbeta(v[k],1,DOF))     
          if(runif(1,0,1)<exp(MHrate)){
            v<-canv;vs<-canvs;probs<-canprobs}
        }
      }
    }    
    
    #update knots:
    for(j in 1:2){for(k in 1:n.terms){
      canvs<-vs;canknot<-knot
      canknot[j,k]<-rnorm(1,knot[j,k],.1)
      if(canknot[j,k]> 0 & canknot[j,k]<1){
        canvs[,k]<-onevs(z,rho,canknot[,k],v[k])
        canprobs<-makeprobs(canvs)
        MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                      log(probs[cbind(1:n,g)]))     
        if(runif(1,0,1)<exp(MHrate)){
          knot<-canknot;vs<-canvs;probs<-canprobs}
      }
    }}    
    
    keepbeta[i,]<-beta
    if(kernel=="gneiting"){
      keepbetagn[i] <- betagn
    } else {
      keeprho[i]<-rho
    }
    
#    count<-count+1
#    if(count==display){
#      par(mfrow=c(2,1))
#      plot(keepbeta[,1])
#      if(kernel=="gneiting"){
#        plot(keepbetagn[1:i])
#      } else {
#        plot(keeprho[1:i])
#      }
#      sdmu<-(mu-min(mu))/(max(mu)-min(mu))
#      plot(z[,1],z[,2],cex=sdmu)
#      plot(knot[,1],knot[,2])
#      count<-0
#    }
    
    if(i>burn){
      afterburn<-afterburn+1
      summu<-summu+mu[g]
      summu2<-summu2+mu[g]^2
      sumtp<-sumtp+probs[,n.terms]
    }
   
    #### Prediction
    if(pred==T){
      vs_new<-matrix(0,nrow(z_new),n.terms)
      if(kernel=="gneiting"){
        for(k in 1:n.terms){vs_new[,k]<-onevs(z_new,rho,knot[,k],v[k],kernel="gneiting",t=t_new,a=1,
                                             beta=betagn)}
      } else {
        for(k in 1:n.terms){vs_new[,k]<-onevs(z=z_new,rho=rho ,
                                              knot = knot[,k],v[k])}
      }
      
      probs_new<-makeprobs(vs_new)
      ypred <- c()
      g_new <- c()
      for(s in 1:nrow(z_new)){
        g_new[s] <- rcat(probs_new[s,])
        ypred[s] <- rnorm(1,x_new[s]*beta + mu[g_new[s]],sige[g_new[s]])
      }
      y_new[i,] <- ypred
    }
  }
  
  post.mn<-summu/afterburn
  post.sd<-summu2/afterburn-post.mn^2
  truncprob<-sumtp/afterburn
  
  list(truncprob=truncprob,beta=keepbeta[burn:runs,],
       probs_last=probs,
       rho=keeprho[burn:runs],betagn=keepbetagn[burn:runs],post.mn=post.mn,post.sd=post.sd,mu_post=mu_post,
       sige_post=sige_post, knot=knot, pred = y_new)}
# 
# 
# load("Y.RData")
# 
# 
# dim(Yout)
# Yout[,,1]
# Yout[,1,]
# Yout[,1,][,5]
# 
# Ydim <- dim(Yout)
# locs <- Ydim[1]
# datcols <- Ydim[3]
# ntimes <- Ydim[2]
# 
# uyears <- unique(Yout[1,,1])
# TT <- length(uyears)
# july <- which((1:ntimes)%%12%%7==0 & (1:ntimes)%%12!=0)
# Y <- t(Yout[,july,5]) # Y is 20 x 100
# 
# normalize <- function(x) {
#   return ((x - min(x)) / (max(x) - min(x)))
# }
# 
# coords <- apply(Yout[,1,][,3:4],2,normalize)
# dim(coords)
# 
#         obj_list_SB <- list()
#         for(j in 1:20){
#           obj_SB <- SSB(y=Y[j,],x=NA,z=coords,
#                       DOF=1,mx.sige=1,mx.sigs=1,n.terms=50,
#                       runs=10000,burn=1000,display=10)
#           obj_list_SB[[uyears[j]]] <- obj_SB
#           save.image("spatialSB.RData")
#         }
# 
# 
# str(obj_list_SB[[1985]])
# clust <- c()
# for(i in 1:100){
#   clust[i] <- which(obj_list_SB[[1985]]$probs_last[i,]==
#                       max(obj_list_SB[[1985]]$probs_last[i,]))
# }
# 
# ylatlon <- Yout[,1,3:4]
# bks=c(14,40)
# quilt.plot(ylatlon[,2],ylatlon[,1],clust,main=uyears[1],
#            ylim=range(ylatlon[,1])+c(-1,1),
#            xlim=range(ylatlon[,2])+c(-1,1))
# map('county',add=T,col='grey')
# map('state',add=T,col='grey60',lwd=2)
# 
#   
#   png("spatialSB_data_85-94")
#   par(mfrow=c(5,2),mar=c(2,3,1,1),fg='grey30',bty='l')
#   for(j in 1:10){
#     clust <- c()
#     for(i in 1:100){
#       clust[i] <- which(obj_list_SB[[uyears[j]]]$probs_last[i,]==
#                         max(obj_list_SB[[uyears[j]]]$probs_last[i,]))
#     }
#     quilt.plot(ylatlon[,2],ylatlon[,1],clust,main=uyears[j],
#                add.legend=F,nlevel=7,nx=45,ny=45,
#                ylim=range(ylatlon[,1])+c(-1,1),
#                xlim=range(ylatlon[,2])+c(-1,1))
#     map('county',add=T,col='grey')
#     map('state',add=T,col='grey60',lwd=2)
#   }
#   dev.off()
#   
#   png("spatialSB_data_95-04")
#   par(mfrow=c(5,2),mar=c(2,3,1,1),fg='grey30',bty='l')
#   for(j in 11:20){
#     clust <- c()
#     for(i in 1:100){
#       clust[i] <- which(obj_list_SB[[uyears[j]]]$probs_last[i,]==
#                           max(obj_list_SB[[uyears[j]]]$probs_last[i,]))
#     }
#     quilt.plot(ylatlon[,2],ylatlon[,1],clust,main=uyears[j],
#                add.legend=F,nlevel=7,nx=45,ny=45,
#                ylim=range(ylatlon[,1])+c(-1,1),
#                xlim=range(ylatlon[,2])+c(-1,1))
#     map('county',add=T,col='grey')
#     map('state',add=T,col='grey60',lwd=2)
#   }
#   dev.off()
#   
#   save.image("spatialSB.RData")
#   