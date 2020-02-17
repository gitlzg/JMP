

library(survival)

#
# Updates the baseline hazard steps using the estimated parameter values
#
updateBaselineHazardSteps = function(censored, surv.times, treatment.status, surv.basis=NA, covariates,
    alphas, phis) {

  options(warn=-1)
  surv.basis=as.matrix(surv.basis)
  nobasis=ncol(surv.basis)
  ddeath=1-as.vector(censored)

  ############################### calculate the unique event times

  data1=cbind(ddeath,surv.times)
  dead.data1=subset(data1,data1[,1]==1)
  uniq.event.times=unique(sort(dead.data1[,2]))
  M=length(uniq.event.times)

  if (!is.na(surv.basis)[1]) {
    #message("dim of surv.basis: ", paste(dim(surv.basis), collapse=","))
    #message("dim of covariates: ", paste(dim(covariates), collapse=","))
    #message("length alphas and phis: ", length(c(alphas, phis)))
    data2=cbind(ddeath,surv.times,surv.basis)
    dead.data2=subset(data2,data2[,1]==1)
    order.surv2 <- dead.data2[order(dead.data2[,2]),]
    uniq.data2=order.surv2[!duplicated(order.surv2[,2]),]
    noba=nobasis+2
    xbeta1=as.vector(treatment.status)%*%(t(uniq.data2[,3:noba,drop=F]%*%alphas))
  } else {
    #message("surv.basis is NA, using just treatment status for baseline hazard steps")
    xbeta1=as.vector(treatment.status)*alphas
  }
  options(warn=0)

  ############################### start to calculate Breslow baseline hazard jumps

  breslow.lam=matrix(NA,M,2)
  breslow.lam[,1]=uniq.event.times

  if (length(phis) > 0) {
    xbeta2=as.vector(covariates%*%phis)
  } else {
    xbeta2 =0
  }

  xbeta=xbeta1+xbeta2
  nob=length(surv.times)
  uniq.times.m=matrix(rep(uniq.event.times,nob),nrow=nob,byrow=T)


  exp.matrix=(as.vector(surv.times)>=uniq.times.m)*exp(xbeta)

  di=(as.vector(surv.times)==uniq.times.m)*ddeath
  dii=apply(di, 2, sum)
  exp.xbeta=apply(exp.matrix, 2, sum)

  breslow.lam[,2]=dii/exp.xbeta

  return(breslow.lam[,2])
}

#
# Function to compute the hazard for specific subject and specific event time.
#
computeHazard = function(baseline.hazard, exp.val) {
  hazard = baseline.hazard * exp.val
  return (hazard)
}

#
# Compute the exponential portion of the hazard. Potentially use the spline basis.
#
computeHazardExp = function(treatment.status, alphas, phis, covariates, surv.basis=NA) {
  exp.val = 0
  betas = c(treatment.status*alphas,phis)
  if (!is.na(surv.basis)[1]) {
    if (!is.matrix(surv.basis)) {
      exp.val = t(c(surv.basis, covariates)) %*% betas
    } else {
      exp.val = cbind(surv.basis, matrix(covariates, nrow=nrow(surv.basis), ncol=length(covariates), byrow=T)) %*% betas
    }
  } else {
    exp.val = t(c(1,covariates)) %*% betas
#    message("Predictors: ", paste(c(1,covariates), collaspe=","))
#    message("betas: ", paste(betas, collaspe=","))
#    message("exp.val ", exp.val)
  }
  return(exp.val)
}

#
# Function to compute the cumulative hazard
#
computeCumulativeHazard = function(surv.time, hazard.exp.val, event.times, baseline.hazard.steps, prior.cum.hazard=NA) {
  cum.hazard = 0
  smaller.uncensored.times = which(event.times <= surv.time)
  largest.smaller.uncensored.time = smaller.uncensored.times[length(smaller.uncensored.times)]
  first.step=1
  if (!is.na(prior.cum.hazard)) {
    first.step = largest.smaller.uncensored.time
    cum.hazard = prior.cum.hazard
  }
  cum.hazard = sum(baseline.hazard.steps[first.step:largest.smaller.uncensored.time]*hazard.exp.val) + cum.hazard
  return (cum.hazard)
}


