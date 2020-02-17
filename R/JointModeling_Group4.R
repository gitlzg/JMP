

# Function to compute the log-likelihood at the specified parameter values for
# all censored subjects without QoL data. This is given by:
#
# -\sum_{i in group4} H(C_i)
#
group4LogLikelihood = function(covariates, # covariate values for just group 4 members
    censored, # indicators of whether each survival time is censored
    all.surv.times, # all survival times (unsorted)
    surv.times, # survival times for just group 4 members
    surv.basis,
    treatment.status, # treatment status for just group 4 members
    alphas, # alpha values (coef. for treatment status)
    phis, # phi values (coef for covariates)
    lambdas # baseline hazard values
) {
  # get unique survival times and sort
  uncensored.event.times = sort(unique(all.surv.times[!censored]))
  if (!is.na(surv.basis)[1]) {
    # subset and reorder the surv.basis for just unique uncensored event times
    indices.to.keep = unique(sapply(uncensored.event.times, function(x) {which(all.surv.times == x)[1]}))
    surv.basis = surv.basis[indices.to.keep,,drop=F]
  }

  logl = 0
  if (!is.matrix(covariates)) { # only one observation so covariates provided as a vector
    covariates = t(covariates)
  }
  n = nrow(covariates) # get the number of subjects in group 4
  # for each subject, compute log-likelihood
  for (i in 1:n) {

    # get the treatment status, survival time, covariate values and baseline hazard at survival for the subject
    treatment.status.i = treatment.status[i]
    surv.time = surv.times[i]
    covariates.i = covariates[i,]
    surv.basis.i = NA
    if(uncensored.event.times[1]<=surv.time){
      if (!is.na(surv.basis)[1]) {
        smaller.uncensored.times = which(uncensored.event.times <= surv.time)
        largest.smaller.uncensored.time = smaller.uncensored.times[length(smaller.uncensored.times)]
        hazard.exp.val = computeHazardExp(treatment.status=treatment.status.i,
            alphas=alphas,
            phis=phis,
            surv.basis=surv.basis[1:largest.smaller.uncensored.time,,drop=F],
            covariates=covariates.i)
      } else {
        hazard.exp.val = computeHazardExp(treatment.status=treatment.status.i,
            alphas=alphas,
            phis=phis,
            surv.basis=NA,
            covariates=covariates.i)
      }
      H_i = computeCumulativeHazard(surv.time=surv.time,
          hazard.exp.val=exp(hazard.exp.val),
          event.times=uncensored.event.times,
          baseline.hazard.steps=lambdas)
    } else {
      # censoring time is before all death times
      H_i = 0
    }

    # Subtract the negative cumulative hazard
    logl = logl - H_i

  }
  return (logl)
}
