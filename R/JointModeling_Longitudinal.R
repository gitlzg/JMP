

#
# Creates the correlation matrix R
#
createR = function(t_i, correlation.exponent) {#}, eta) {
  n_i = length(t_i)
  return (diag(1, n_i))
#  R = matrix(0, nrow=n_i, ncol=n_i)
#  for (j in 1:n_i) {
#    for (l in 1:n_i) {
#      R[j,l] = exp(-eta*abs(t_i[j] -t_i[l])^correlation.exponent)
#    }
#  }
#  return (R)
}

#
# Compute the longitudinal mean a specific time and subject
# given the specified mu and beta param values
#
computeLongitudinalMean = function(long.basis, mus, betas, treatment.status) {

  # compute the value of the mu spline
  mu.value = cubicSpline(long.basis, params=mus)

  # if treated, compute the value of the beta spline
  if (treatment.status) {
    beta.value = cubicSpline(long.basis, params=betas)
    value = mu.value + beta.value
  } else {
    value = mu.value
  }

  return (as.vector(value))
}

#
# Creates the cache for computation of the R matrix
#
createRCache = function() {
  cache = new.env(parent=globalenv())
  cache$saved_nu = -1
  cache$saved_sigma = -1
  #cache$saved_eta = -1
  cache$saved_R_i = list()
  cache$saved_det.V_i = list()
  cache$saved_inv.V_i = list()
  cache$saved_long_basis = NA
  cache$saved_uncensored.event.times = NA
  cache$saved_surv.time.for.event = NA
  cache$eval.num = 0
  return (cache)
}

#
# JIT creation of R_i related objects
#
populateRCache = function(long.status, num.long.measurements, long.times,
     correlation.exponent, sigma, nu, cache) {
   #eta, cache) {
  #if (sigma != cache$saved_sigma | eta != cache$saved_eta | nu != cache$saved_nu) {
  if (sigma != cache$saved_sigma | nu != cache$saved_nu) {
    #message("Populating R cache...")
    n = nrow(long.status)
    # JIT R_i creation
    #if (eta != cache$saved_eta) {
    # cache$saved_eta = eta
      for (i in 1:n) {
        n_i = num.long.measurements[i]
        long.measurements = which(long.status[i,] ==1)
        t_i = as.vector(long.times[i,long.measurements])
        cache$saved_R_i[[i]] = createR(t_i, correlation.exponent)#, eta)
      }
    #}
    cache$saved_sigma = sigma
    cache$saved_nu = nu
    for (i in 1:n) {
      n_i = num.long.measurements[i]
      R_i = cache$saved_R_i[[i]]
      #message("Length n_i: ", n_i, ", nrow R_i: ", nrow(R_i))
      #message("Computing V_i using sigma: ", sigma)
      #V_i = tau^2*diag(rep(1,n_i)) + sigma^2*matrix(1,nrow=n_i, ncol=n_i) + nu^2*R_i
      V_i =  sigma^2*matrix(1,nrow=n_i, ncol=n_i) + nu^2*R_i
      cond.num = rcond(V_i)
      #message("V_i condition number: ", cond.num)
      if (!is.nan(cond.num) & cond.num >= .Machine$double.eps) {
        #message("About to get compute det/inverse...")
        cache$saved_det.V_i[[i]] = det(V_i)
        cache$saved_inv.V_i[[i]] = solve(V_i)
        #message("...computed det/inverse.")
      } else {
        message("V_i is near singular, condition number: ", cond.num)
        cache$saved_det.V_i[[i]] = NA
        cache$saved_inv.V_i[[i]] = NA
      }
    }
  }
}
