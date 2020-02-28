
#-------------------------------------------------------------------
# Functions for computing the group 1 and 2 survival likelihood
#-------------------------------------------------------------------

#library(survival)

#
# Function to compute the log-likelihood at the specified parameter values for
# all uncensored subjects without QoL data.
#
group12SurvivalLogLikelihood = function(covariates, # covariate values for just group 1 and 2 members
    surv.times, # survival times for just group 1 and 2 members, i.e., just the uncensored times
    surv.basis=NA,
    treatment.status, # treatment status for just group 1 and 2 members
    alphas, # alpha values
    phis, # phi values
    lambdas # baseline hazard values
) {
  # get unique survival times and sort
  unique.surv.times = unique(surv.times)
  unique.event.times = sort(unique.surv.times)

  if (!is.na(surv.basis)[1]) {
    # subset and reorder the surv.basis for just unique uncensored event times
    indices.to.keep = unique(sapply(unique.event.times, function(x) {which(surv.times == x)[1]}))
    surv.basis = surv.basis[indices.to.keep,,drop=FALSE]
  }

  # for each group member, determine which of the unique sorted events corresponds to their survival time
  hazard.step.per.subject = sapply(surv.times, function(x) {which(unique.event.times == x)})
  logl = 0

  if (!is.matrix(covariates)) { # only one observation so covariates provided as a vector
    covariates = t(covariates)
  }

  # determine the number of group 1 and 2 members
  n = nrow(covariates)

  # For each group member, compute the log-likelihood contribution
  for (i in 1:n) {

    # get the treatment status, survival time, covariate values and baseline hazard at survival for the subject
    treatment.status.i = treatment.status[i]
    surv.time = surv.times[i]
    covariates.i = covariates[i,]
    surv.basis.i = NA
    if (!is.na(surv.basis)[1]) {
      surv.basis.i= surv.basis[hazard.step.per.subject[i],]
    }
    baseline.hazard = lambdas[hazard.step.per.subject[i]]

    # compute the hazard for the subject
    hazard.exp.val = computeHazardExp(treatment.status=treatment.status.i,
        alphas=alphas,
        phis=phis,
        covariates=covariates.i,
        surv.basis=surv.basis.i)
    hazard = computeHazard(baseline.hazard=baseline.hazard,
            exp.val=exp(hazard.exp.val))
    #message("hazard for subject ", i, ": ", hazard)

    # get the log of the hazard
    logl_i = log(hazard)

    # if there is a survival spline, need to compute the exponential hazard at each time
    if (!is.na(surv.basis)[1]) {
      smaller.uncensored.times = which(unique.event.times <= surv.time)
      largest.smaller.uncensored.time = smaller.uncensored.times[length(smaller.uncensored.times)]
      hazard.exp.val = computeHazardExp(treatment.status=treatment.status.i,
          alphas=alphas,
          phis=phis,
          surv.basis=surv.basis[1:largest.smaller.uncensored.time,,drop=FALSE],
          covariates=covariates.i)
    }

    # subtract the cumulative hazard
    logl_i = logl_i - computeCumulativeHazard(surv.time=surv.time,
        hazard.exp.val=exp(hazard.exp.val),
        event.times=unique.event.times,
        baseline.hazard.steps=lambdas)
    #message("Survival likelihood for subject ", i, ": ", logl_i)
    logl = logl + logl_i
  }
  return (logl)
}
