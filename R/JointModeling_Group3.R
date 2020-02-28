

#-------------------------------------------------------------------
# Functions for computing the group 3 likelihood
#-------------------------------------------------------------------

populateBasisCache = function(n, surv.times, censored, long.status, long.times, long.knots, surv.times.for.group, cache,
    natural.spline=FALSE) {
  options(warn=-1)
  if (!is.na(cache$saved_uncensored.event.times)) {
    #message("basis cache already populated")
    return ()
  } else {
    #message("need to populate cache")
  }
  options(warn=0)

  cache$saved_uncensored.event.times = sort(unique(surv.times[!censored]))
  cache$saved_surv.time.for.event = sapply(cache$saved_uncensored.event.times, function(x) {which(surv.times == x)[1]})
  cache$saved_long.basis = list()
  for (i in 1:n) {
    cache$saved_long.basis[[i]] = list()
    surv.time = surv.times.for.group[i]
    long.measurements = which(long.status[i,] ==1)
    t_i = long.times[i,long.measurements]
    for (evt in 1:length(cache$saved_uncensored.event.times)) {
      evt.time = cache$saved_uncensored.event.times[evt]
      if (evt.time <= surv.time) {
        # only process event times greater than survival time
        next
      }
      #message("Event time: ", evt.time, ", cenored time: ", surv.time, ", measurement times: ", paste(t_i, collapse=", "))
      adj.times = sapply(t_i, function(x) {
            evt.time - (surv.time - x)
          })
      cache$saved_long.basis[[i]][[evt]] = getSplineBasis(knots=long.knots, times=adj.times, natural.spline=natural.spline)
    }
  }
  #message("Finished with populationBasisCache()")
}

#
# Function to compute the log-likelihood at the specified parameter values for
# all censored subjects with QoL data.
#
group3LogLikelihood = function(long.status, num.long.measurements,
    long.values, long.times,
    long.covariates, surv.covariates,
    treatment.status, long.knots, natural.spline,
    surv.basis,
    correlation.exponent,
    mus, betas, psis, sigma, nu, cache, #eta, cache,
    surv.times, surv.times.for.group, censored,
    alphas, phis, lambdas,
    special.weighting=TRUE,
    most.likely.event=FALSE,
    long.log.sum.exp.trick=TRUE,
    haz.log.sum.exp.trick=FALSE) {

  logl = 0
  if (!is.matrix(long.covariates)) { # only one observation so everything is provided as a vector
    long.covariates = t(long.covariates)
    surv.covariates = t(surv.covariates)
    long.status = t(long.status)
    long.values = t(long.values)
    long.times = t(long.times)
  }

  n = nrow(long.covariates) # get the number of subjects in group 3

  populateRCache(long.status=long.status,
      num.long.measurements=num.long.measurements,
      long.times=long.times,
      correlation.exponent=correlation.exponent,
      sigma=sigma,
      nu=nu,
      #eta=eta,
      cache=cache)

  populateBasisCache(n, surv.times, censored,
      long.status, long.times, long.knots, surv.times.for.group,
      cache, natural.spline)

  if (ncol(long.covariates) > 0) {
    cov_psis = long.covariates %*% psis
  } else {
    cov_psis = rep(0,n)
  }

  if (!is.na(surv.basis)[1]) {
    # subset and reorder the surv.basis for just unique uncensored event times
    unique.surv.times = unique(surv.times)
    indices.to.keep = unique(sapply(cache$saved_uncensored.event.times, function(x) {which(surv.times == x)[1]}))
    surv.basis = surv.basis[indices.to.keep,,drop=FALSE]
  }

  # Compute likelihood for each group 3 subject
  for (i in 1:n) {

    #message("Computing likelihood for group 3 subject ", i, "...")
    treatment.status.i = treatment.status[i]
    surv.time = surv.times.for.group[i]
    long.covariates.i = long.covariates[i,]
    surv.covariates.i = surv.covariates[i,]
    n_i = num.long.measurements[i]
    long.measurements = which(long.status[i,] ==1)
    t_i = long.times[i,long.measurements]
    Y_i = long.values[i,long.measurements]
    R_i = cache$saved_R_i[[i]]
    det.V_i = cache$saved_det.V_i[[i]]
    if (is.na(det.V_i)) {  # V_i was singular
      return (NaN)
    }
    inv.V_i = cache$saved_inv.V_i[[i]]
    # for each subject, sum over all event times that are >= the censoring time
    l_for_subject = 0
    cum.hazard = NA
    if (!is.na(surv.basis)[1]) {
      #alphas = rep(0, length(alphas))
      all.hazard.exp.val = computeHazardExp(treatment.status=treatment.status.i,
          alphas=alphas,
          phis=phis,
          surv.basis=surv.basis,
          covariates=surv.covariates.i)
    }

    #message("Computing cumulative hazard, total events: ", length(event.times), ", subject event: ", surv.hazard.step)

    # if special weights are being used or only the most likely future death time or the log-sum-exp trick for the survival portion of the
    # likelihood, need to precompute various items
    haz.exp.cum.haz.values = NA
    most.likely.event.index = 0
    max.neg.cum.hazard = 0
    S_i = NA

    # if using special weights, just the most likely event of the log-sum-exp trick for the survival portion
    # of the likelihood, need to precompute various items
    if (special.weighting | most.likely.event | haz.log.sum.exp.trick) {
      S_i = 0
      cum.hazard = 0
      censoring.cum.hazard = 0
      cum.hazards = rep(NA, length(cache$saved_uncensored.event.times))
      later.cum.hazards = rep(NA, length(cache$saved_uncensored.event.times))
      haz.exp.cum.haz.values = rep(0, length(cache$saved_uncensored.event.times))

      # compute the cumulative hazard at all death times; save the cum hazard at the
      # censoring time
      for (evt in 1:length(cache$saved_uncensored.event.times)) {
        evt.time = cache$saved_uncensored.event.times[evt]
        cum.hazard = cum.hazard +
            computeHazard(baseline.hazard=lambdas[evt], exp.val=exp(all.hazard.exp.val[evt]))
        cum.hazards[evt] = cum.hazard
        evt.time = cache$saved_uncensored.event.times[evt]
        if (evt.time < surv.time) {
          censoring.cum.hazard = cum.hazard
        } else if (evt.time > surv.time) {
          later.cum.hazards[evt] = cum.hazard
        }
      }

      # For the log-sum-exp trick, need the max of the hazard exp portion and
      # the negative cum hazard for all death times later than the censoring time
      max.hazard.exp.val = max(all.hazard.exp.val[which(!is.na(later.cum.hazards))])
      max.neg.cum.hazard = max(-later.cum.hazards[which(!is.na(later.cum.hazards))])

      # Compute and cache the survival density at all death times after the censoring time
      # Also calculate the weight, S_i
      #message("Subject ", i, ", largest event time: ", cache$saved_uncensored.event.times[length(cache$saved_uncensored.event.times)],
      #    ", surv time: ", surv.time)
      for (evt in 1:length(cache$saved_uncensored.event.times)) {
        evt.time = cache$saved_uncensored.event.times[evt]
        if (evt.time <= surv.time) {
          # only process event times greater than survival time
          next
        }
        if (haz.log.sum.exp.trick) {
          hazard = computeHazard(baseline.hazard=lambdas[evt], exp.val=exp(all.hazard.exp.val[evt]-max.hazard.exp.val))
          haz.exp.cum.haz.values[evt] = hazard*exp(-cum.hazards[evt] - max.neg.cum.hazard)
        } else {

          hazard = computeHazard(baseline.hazard=lambdas[evt], exp.val=exp(all.hazard.exp.val[evt]))

          haz.exp.cum.haz.values[evt] = hazard*exp(-cum.hazards[evt])
          #message("Adding to S_i for subject ", i, ": ", haz.exp.cum.haz.values[evt], ", hazard: ", hazard, ", exp(-cum.hazard): ", exp(-cum.hazards[evt]))
        }

        S_i = S_i + haz.exp.cum.haz.values[evt]
      }
      #message("S_i: ", S_i, ", exp(-H_c): ", exp(-censoring.cum.hazard))
      S_i = S_i/exp(-censoring.cum.hazard)

      # Find the event with the largest density
      if (most.likely.event) {
        next.largest.event = 0
        for (evt in 1:length(cache$saved_uncensored.event.times)) {
          evt.time = cache$saved_uncensored.event.times[evt]
          if (evt.time <= surv.time) {
            # only process event times greater than survival time
            next
          }
          next.largest.event = evt
          break
        }
        larger.haz.exp.cum.haz.values = haz.exp.cum.haz.values[next.largest.event:length(cache$saved_uncensored.event.times)]
        most.likely.event.index = which(larger.haz.exp.cum.haz.values == max(larger.haz.exp.cum.haz.values)) + next.largest.event -1
       # message("Most likely future death event index is ", most.likely.event.index,
       #     " from events ", next.largest.event, " to ", length(cache$saved_uncensored.event.times))
      }
    }

    # Evaluate likelihood at each future death time (potentially at just the
    # most likely). This is looped through twice to support the
    # log-sum-exp trick for longitudinal portion of the likelihood

    total.evaluated.events = 0
    #long.recip.sqrt.det = 1/sqrt(det.V_i)
    long.exps = rep(NA, length(cache$saved_uncensored.event.times))

    for (evt in 1:length(cache$saved_uncensored.event.times)) {
      evt.time = cache$saved_uncensored.event.times[evt]
      if (evt.time <= surv.time) {
        # only process event times greater than survival time
        next
      }
      total.evaluated.events = total.evaluated.events +1
      if (most.likely.event) {
        if (evt != most.likely.event.index) {
          # only process most likey next death time
          next
        }
      }
      if (i > length(cache$saved_long.basis)) {
        message("i: ", i, ", length saved_long.basis: ", length(cache$saved_long.basis))
      } else if (evt > length(cache$saved_long.basis[[i]])) {
        message("evt: ", evt, ", length saved_long.basis[[i]]: ",
            length(cache$saved_long.basis[[i]]))
      }

      long.basis = cache$saved_long.basis[[i]][[evt]]
      f_i = computeLongitudinalMean(long.basis, mus, betas, treatment.status.i)
      long.delta = Y_i - f_i - cov_psis[i]
      long.exps[evt] = -.5*(t(long.delta) %*% inv.V_i %*% long.delta)
      #long.recip.sqrt.det = 1/((2*pi)^(n_i/2)*sqrt(det.V_i))
    }

    a = 0
    # If using the log-sum-exp trick, compute the constant
    if (long.log.sum.exp.trick) {
      a = max(long.exps[which(!is.na(long.exps))])
    }

    # debugging structure
    long.densities = rep(0, length(cache$saved_uncensored.event.times))

    for (evt in 1:length(cache$saved_uncensored.event.times)) {
      evt.time = cache$saved_uncensored.event.times[evt]
      if (evt.time <= surv.time) {
        # only process event times greater than survival time
        next
      }
      if (most.likely.event) {
        if (evt != most.likely.event.index) {
          # only process most likey next death time
          next
        }
      }
      long.exp = long.exps[evt]
      # Using the log-sum-exp trick: https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/
      if (long.log.sum.exp.trick) {
        #long.density = long.recip.sqrt.det*exp(long.exp-a)
        long.density = exp(long.exp-a)
      } else {
        #long.density = long.recip.sqrt.det*exp(long.exp)
        long.density = exp(long.exp)
      }
      long.densities[evt] = long.density

      if (special.weighting | most.likely.event| haz.log.sum.exp.trick) {
        # if special weighting is being used, these have already been computed
        haz.exp.cum.haz = haz.exp.cum.haz.values[evt]
      } else {
        cum.hazard = computeCumulativeHazard(surv.time=evt.time,
            hazard.exp.val=exp(all.hazard.exp.val[1:evt]),
            event.times=cache$saved_uncensored.event.times[1:evt],
            baseline.hazard.steps=lambdas[1:evt])
        hazard = computeHazard(baseline.hazard=lambdas[evt], exp.val=exp(all.hazard.exp.val[evt]))
        haz.exp.cum.haz = hazard*exp(-cum.hazard)
      }

      # add on the likelihood for this subject for this event time
      if (special.weighting) {
        if (haz.exp.cum.haz == 0) { # protect against 0/0
          evt.density = 0
        } else {
          evt.density = long.density*(haz.exp.cum.haz/S_i)
        }
      } else {
        evt.density = long.density*haz.exp.cum.haz
      }

      #message("Density for event time ", evt.time, " for subject ", i, " with surv time: ", surv.time, ": ", evt.density,
      #  ", long.density: ", long.density, ", hazard: ", hazard, ", cum.hazard: ", cum.hazard,
      #  ", det.V_i: ", det.V_i, ", exp: ", long.exp, ", long.delta: ", paste(long.delta, collapse=","))
      l_for_subject = l_for_subject + evt.density
      # IMPORTANT: uncomment this break to just compute the likelihood at the first death time after
      #            censoring time; this will make group 3 much faster for testing
      #break
    }

    # if the likelihood for the subject is below threshold even after log-sum-exp trick, output the associated values
    if (is.infinite(log(l_for_subject))) {
      message("Likelihood for group 3 subject ", i, " is below minimum: ", l_for_subject,
          ", long.exp: ", long.exp, ", long.density: ", long.density,
          ", a: ", a,
          ", long.exps: ", paste(long.exps, collapse=","),
          ", length(long.exps): ", length(long.exps),
          ", haz.exp.cum.haz: ", haz.exp.cum.haz,
          ", S_i: ", S_i,
          ", long.delta: ", (t(long.delta)%*% long.delta),
          ", haz.exp.cum.haz: ", haz.exp.cum.haz,
          ", cum.hazards: ", paste(cum.hazards, collapse=", "),
          ", haz.exp.cum.haz.values: ", paste(haz.exp.cum.haz.values, collapse=","),
          ", long.densities: ", paste(long.densities, collapse=","))
    }
    # take the log
    logl_for_subject = log(l_for_subject)-0.5*log(det.V_i)
    if (long.log.sum.exp.trick) {
      logl_for_subject = logl_for_subject+a
    }
    if (haz.log.sum.exp.trick & !special.weighting) {
        logl_for_subject = logl_for_subject+max.neg.cum.hazard+max.hazard.exp.val
    }

    #message("Computed likelihood for group 3 subject ", i, ": ", logl_for_subject)
    # add on the likelihood for the subject
    logl = logl + logl_for_subject
  }

  #message("Computed group 3 likelihood: ", logl)

  return (logl)
}
