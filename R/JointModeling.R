

library(numDeriv)



#-------------------------------------------------------------------
# Main public function to maximize the joint longitudinal and survival
# likelihood. Can also be used to just compute the log-likelihood at a specific set of
# parameter values.
#-------------------------------------------------------------------
optimizeJointLikelihood = function(
    long.values, # matrix of longitudinal values, one row per subject, one column per measurement
    long.times, # matrix of longitudinal measurement times in months from death
    surv.times, # vector of survival time (potentially censored) in months from enrollment
    censored, # vector of booleans indicating censoring status
    treatment.status, # vector of binary treatment status (1 for treated)
    long.covariates=NA, # matrix of fixed effect covariates, NA for no covariates
    surv.covariates=NA, # matrix of fixed effect covariates, NA for no covariates
    long.knots=c(), # vector of cubic spline knot times for the longitudinal model (from death in months)
    surv.knots=c(), # vector of cubic spline knot times for the survival model (from enrollment in months)
    natural.spline=FALSE,
    correlation.exponent=2, # the correlation exponent for the R matrix (1 or 2)
    param.init.values, # initial values of the parameters, MUST be specified in the following order
    # mu (mean longitudinal value, cubic spline basis params)
    # beta (coef of treatment indicator in longitudinal model, cubic spline basis params)
    # psi (coef of covariates in longitudinal model, vector of params)
    # sigma (random effect covariance matrix)
    # nu (stddev of longitudinal Gaussian process)
    # NOT INCLUDED eta (exponential coef in Gaussian process correlation)
    # alpha (coef of treatment indicator in survival model)
    # phi (coef of covariates in survival model, vector)
    param.lower.bounds, # lower bounds on the parameters
    param.upper.bounds, # upper bounds on the parameters
    optim.method=NA, # optimization method, NA to just compute the log-likelihood at the initial parameter values
    control=NA, # control object for calling optim()
    outer.iterations=100, # number of outer iterations for calls to constrOptim(); done for all
    # methods except "L-BFGS-B"
    outer.eps=1e-5, # convergence tolerance for constrOptim()
    enforce.bounds=TRUE, # True to enforce upper and lower param bounds via constrOptim(), false to just call optim()
    include.survival=TRUE, # boolean controling whether the likelihood for the survival model is computed
    include.longitudinal= TRUE, # boolean controling whether the likelihood for the longitudinal model is computed
    compute.standard.errors=FALSE, # True to compute the Hessian
    profile=FALSE, # True to profile the execution time
    naive=FALSE # True to treat group 3 as group 4
  ){

  # Compute some needed items
  num.events = length(unique(surv.times))
  num.long.covariates = 0
  if (length(long.covariates) > 0) {
    num.long.covariates = ncol(long.covariates)
  }
  num.surv.covariates = 0
  if (length(surv.covariates) > 0) {
    num.surv.covariates = ncol(surv.covariates)
  }

  hess = NA

  param.indexes = getParamIndexes(
      num.long.fixed.effects=num.long.covariates,
      num.surv.fixed.effects=num.surv.covariates,
      num.long.knots=length(long.knots),
      num.surv.knots=length(surv.knots),
      num.events=num.events,
      include.survival=include.survival,
      include.longitudinal=include.longitudinal,
      natural.spline=natural.spline)

  if (include.longitudinal) {
    long.status = apply(long.values, c(1,2), function(x) {as.numeric(!is.na(x))})
    long.time.status = apply(long.times, c(1,2), function(x) {as.numeric(!is.na(x))})
    long.status = long.status + long.time.status
    long.status = apply(long.status, c(1,2), function(x) {if (x == 2) {return (1)} else {return(0)}})
    num.long.measurements = apply(long.status, 1, sum)
  } else {
    long.status = NA
    num.long.measurements = 0
  }

  # Determine membership of subjects in the 4 groups:
  has.long.data = num.long.measurements > 0
  group1 = which(!censored & has.long.data)
  group2 = which(!censored & !has.long.data)
  group3 = which(censored & has.long.data)
  group4 = which(censored & !has.long.data)

  # get all unique longitudinal measurement times for group 1
  unique.times = unique(as.vector(as.matrix(long.times[group1,])))

  # remove any NA values
  unique.times = unique.times[which(!is.na(unique.times))]

  # create a matrix of the same size as long times but with the indexes
  # of the time in the vector of unique times
  if (include.longitudinal) {
    unique.time.indexes = matrix(0, nrow=length(group1), ncol=ncol(long.times))
    for (i in 1:length(group1)) {
      long.measurements = which(long.status[group1,][i,] ==1)
      if (length(long.measurements) > 0) {
        times = as.vector(long.times[group1,][i,long.measurements])
        indexes = sapply(times, function(x){which(unique.times == x)})
        unique.time.indexes[i,long.measurements] = indexes
      }
    }
    long.basis = getSplineBasis(knots=long.knots, times=unique.times, natural.spline=natural.spline)
    message("Num cols long.basis: ", ncol(long.basis), ", num betas: ", length(param.indexes$beta.indexes))
  } else {
    unique.time.indexes=NA
    long.basis=NA
  }

  # create the survival basis, get it for all subject surv.times in the original order
  if (include.survival) {
    surv.basis = getSplineBasis(knots=surv.knots, times=surv.times, natural.spline=natural.spline)
    message("Num cols surv.basis: ", ncol(surv.basis), ", num alphas: ", length(param.indexes$alpha.indexes))
  } else {
    surv.basis=NA
  }

  # Create a list of arguments for the joint log-likelihood function (for code cleanliness)
  args = list(long.status=long.status,
    num.long.measurements=num.long.measurements,
    long.values=long.values,
    long.times=long.times,
    long.covariates=long.covariates,
    surv.covariates=surv.covariates,
    surv.times=surv.times,
    treatment.status=treatment.status,
    censored=censored,
    long.knots=long.knots,
    surv.knots=surv.knots,
    natural.spline=natural.spline,
    correlation.exponent=correlation.exponent,
    mu.indexes=param.indexes$mu.indexes,
    beta.indexes=param.indexes$beta.indexes,
    psi.indexes=param.indexes$psi.indexes,
    sigma.index=param.indexes$sigma.index,
    tau.index=param.indexes$tau.index,
    nu.index=param.indexes$nu.index,
    #eta.index=param.indexes$eta.index,
    alpha.indexes=param.indexes$alpha.indexes,
    phi.indexes=param.indexes$phi.indexes,
    include.survival=include.survival,
    include.longitudinal=include.longitudinal,
    long.basis=long.basis,
    surv.basis=surv.basis,
    unique.time.indexes=unique.time.indexes,
    iteration=1,
    naive=naive)

  # Cache of saved data structures
  group1.R.cache = createRCache()
  group3.R.cache = createRCache()
  if (is.na(optim.method)) {

    message("Computing log-likelihood without optimization...")

    # try generating the likelihood once
    num.executions = 1
    if (profile) {
      num.executions=10
      Rprof("JointModelingProf.out")
    }
    start = proc.time()[3]
    for (i in 1:num.executions) {
    neg.logl = jointLogLikelihood(params=param.init.values, args=args,
        group1.R.cache=group1.R.cache,
        group3.R.cache=group3.R.cache)
    }
    execution.time = (proc.time()[3] - start)/10
    if (profile) {
      Rprof(NULL)
    }

    # optionally compute Hessian and standard errors
    if (compute.standard.errors) {
      hess = jointLogLikelihoodHessian(params=param.init.values,
          args=args,
          group1.R.cache=group1.R.cache,
          group3.R.cache=group3.R.cache)
    }

    # build the result list
    result = list()
    result$param.init.values = param.init.values
    result$param.lower.bounds=param.lower.bounds
    result$param.upper.bounds=param.upper.bounds
    result$neg.log.likelihood = neg.logl
    result$hessian=hess
    result$standard.errors = NA
    result$execution.time = execution.time
    if (!is.na(hess)) {
      result$standard.errors = sqrt(diag(solve(hess)))
    }

    return (result)

  } else {

    message("Optimizing log-likelihood...")

    # Call constrOptim
    constraints = buildConstOptimStructures(param.upper.bounds, param.lower.bounds)

    # check initial value against constraints
    delta = constraints$ui %*% param.init.values - constraints$ci
    violated.constraints = which(delta < 0)
    if (length(violated.constraints) > 0) {
      for (i in 1:length(violated.constraints)) {
        violated.constraint = violated.constraints[i]
        message("Violated constraint ", violated.constraint, ", ui: ", paste(constraints$ui[violated.constraint,], collapse=" "),
            ", ci: ", constraints$ci[violated.constraint], ", delta: ", delta[violated.constraint])
      }
      stop("Constraints were violated by initial parameter values!")
    }

    if (optim.method=="SANN") {
      grad.method = NULL
    } else {
      grad.method = jointLogLikelihoodGrad
    }

    start = proc.time()[3]
    args$iteration="Optimize and Hessian"

    if (enforce.bounds) {
      optim.results = constrOptim(theta=param.init.values,
         f=jointLogLikelihood,
         grad=grad.method,
         ui=constraints$ui,
         ci=constraints$ci,
         method=optim.method,
         outer.iterations=outer.iterations,
         outer.eps=outer.eps,
         control=control,
         args=args,
         group1.R.cache=group1.R.cache,
         group3.R.cache=group3.R.cache)
     if (compute.standard.errors) {
       message("Computing joint log-likelihood hessian")
       optim.results$hessian = jointLogLikelihoodHessian(params=optim.results$par,
           args=args, group1.R.cache=group1.R.cache, group3.R.cache=group3.R.cache)
     }
    } else {
     grad = grad.method(param.init.values, args, group1.R.cache,group3.R.cache)
     message("sum of gradients: ", sum(grad))
     ll = jointLogLikelihood(param.init.values, args, group1.R.cache,group3.R.cache)
     message("Likelihood: ", ll)
     message("About to call optim with initial values: ", paste(param.init.values, collapse=","))
     optim.results = optim(par=param.init.values,
         fn=jointLogLikelihood,
         gr=grad.method,
         args=args,
         group1.R.cache=group1.R.cache,
         group3.R.cache=group3.R.cache,
         method=optim.method,
         control=control)
     if (compute.standard.errors) {
       optim.results$hessian = jointLogLikelihoodHessian(params=optim.results$par,
           args=args, group1.R.cache=group1.R.cache, group3.R.cache=group3.R.cache)
     }
     message("Called optim")
    }

    execution.time = proc.time()[3] - start
    if (compute.standard.errors) {
      hess = optim.results$hessian
    }

    result = list()
    result$param.init.values = param.init.values
    result$param.lower.bounds=param.lower.bounds
    result$param.upper.bounds=param.upper.bounds
    result$optim.results = optim.results
    result$standard.errors = NA
    # compute the standard errors from the hessian (since negative hessian is returned
    # when minimizing -log-likelihood, can directly invert)
    if (!is.na(hess)) {
      result$cov.mat = solve(hess)
      result$standard.errors = sqrt(diag(result$cov.mat))
    }
    result$optim.method = optim.method
    result$neg.log.likelihood = optim.results$value
    result$optim.control = control
    result$execution.time = execution.time

    return (result)
   }
 }


#-------------------------------------------------------------------
# Create the u_i matrix and c_i vector for use with constrOptim
#-------------------------------------------------------------------
buildConstOptimStructures = function(upper, lower) {
  p = length(upper)
  ui = matrix(0, nrow=2*p, ncol=p)
  ci = rep(0, nrow=p)
  for (i in 1:p) {
    # specify the min
    ui[1+(i-1)*2,i] = 1
    ci[1+(i-1)*2] = lower[i]
    # specify the max
    ui[2+(i-1)*2,i] = -1
    ci[2+(i-1)*2] = -upper[i]
  }
  return (list(ui=ui, ci=ci))
}

#-------------------------------------------------------------------
# Define the joint likelihood gradient function
#-------------------------------------------------------------------
jointLogLikelihoodGrad = function(params, args,
    group1.R.cache=group1.R.cache,
    group3.R.cache=group3.R.cache) {
  # defaults to Richardson method
  return (grad(func=jointLogLikelihood, params,
          method="simple",
          args=args,
          group1.R.cache=group1.R.cache,
          group3.R.cache=group3.R.cache))
}

#-------------------------------------------------------------------
# Compute the joint likelihood hessian function
#-------------------------------------------------------------------
jointLogLikelihoodHessian = function(params, args,
    group1.R.cache=group1.R.cache,
    group3.R.cache=group3.R.cache) {
  return (jacobian(func=jointLogLikelihoodGrad,
      x=params,
      method="simple",
      args=args,
      group1.R.cache=group1.R.cache,
      group3.R.cache=group3.R.cache))
}



