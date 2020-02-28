
#-------------------------------------------------------------------
# Functions for computing the group 1 longitudinal likelihood
#-------------------------------------------------------------------


#
# Function to compute the log-likelihood at the specified parameter values for
# all uncensored subjects with longitudinal data.
#
group1LongitudinalLogLikelihood = function(long.status, num.long.measurements,
    long.values, long.times,
    covariates, treatment.status, long.basis, unique.time.indexes, correlation.exponent,
    mus, betas, psis, sigma, nu, cache){ #eta, cache) {

  logl = 0

  if (!is.matrix(covariates)) { # only one observation so everything is provided as a vector
    covariates = t(covariates)
    long.status = t(long.status)
    long.values = t(long.values)
    long.times = t(long.times)
  }

  n = nrow(long.status)

  populateRCache(long.status=long.status,
      num.long.measurements=num.long.measurements,
      long.times=long.times,
      correlation.exponent=correlation.exponent,
      sigma=sigma,
      nu=nu,
      #eta=eta,
      cache=cache)

  for (i in 1:n) {
    n_i = num.long.measurements[i]
    long.measurements = which(long.status[i,] ==1)
    time.indexes = as.vector(unique.time.indexes[i, long.measurements])
    Y_i = long.values[i,long.measurements]
    cov_i = covariates[i,]
    R_i = cache$saved_R_i[[i]]
    det.V_i = cache$saved_det.V_i[[i]]
    if (is.na(det.V_i) | is.infinite(det.V_i)) {  # V_i was singular
      return (NaN)
    }
    inv.V_i = cache$saved_inv.V_i[[i]]
    cov_i_psi = 0
    if (length(cov_i) > 0) {
      cov_i_psi = t(cov_i) %*% psis
    }

    subset.long.basis = long.basis[time.indexes,,drop=FALSE]
    f_i = computeLongitudinalMean(subset.long.basis, mus, betas, treatment.status[i])
    long.delta = suppressWarnings(Y_i - f_i - cov_i_psi)
    logl_i = -.5*(n_i*log(2*pi) + log(det.V_i) + t(long.delta) %*% inv.V_i %*% long.delta)
#    if (is.na(logl_i)) {
#      message("NA log-likelihood for group 1 subject ", i, ", skipping. ",
#          "Y_i: ", paste(Y_i, collapse=","),
#          "f_i: ", paste(f_i, collapse=","),
#          "cov_i_psi: ", paste(cov_i_psi, collapse=","),
#          "cov_i: ", paste(cov_i, collapse=","),
#          "psis: ", paste(psis, collapse=","),
#          ", long.delta: ", long.delta, ", det.V_i: ", det.V_i, ", n_i: ", n_i)
#      next
#    }
    logl = logl + logl_i
  }
  return (logl)
}
