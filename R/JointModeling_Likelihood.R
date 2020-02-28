
#------------------------------------------------------------------------------------
# Function to compute the joint log-likelihood of the longitudinal
# and survival models given the empirical data and specified values of the parameters.
# Returns the -log-likelihood so that optimization can be performed via minimization.
#------------------------------------------------------------------------------------
jointLogLikelihood = function(
    params, # vector of parameter values
    args, # list of arguments for computing log-likelihood (built in optimizeJointLikelihood())
    group1.R.cache,
    group3.R.cache
){

  # if (group1.R.cache$eval.num %% 50 == 0) {
  #   message("Evaluation ", group1.R.cache$eval.num, ", computing log-likelihood with params: ", paste(params, collapse=", "))
  # }

  # Recompute the baseline hazard steps using the new parameter values
  if (args$include.survival) {
    #message("Updating baseline hazard")
    lambdas = updateBaselineHazardSteps(censored=args$censored,
        surv.times=args$surv.times,
        treatment.status=args$treatment.status,
        surv.basis=args$surv.basis,
        covariates=args$surv.covariates,
        alphas=params[args$alpha.indexes],
        phis=params[args$phi.indexes])
  }

  # Determine membership of subjects in the 4 groups:
  has.long.data = args$num.long.measurements > 0
  group1 = which(!args$censored & has.long.data)
  group2 = which(!args$censored & !has.long.data)
  group3 = which(args$censored & has.long.data)
  if (args$naive) {
    group4 = which(args$censored)
  } else {
    group4 = which(args$censored & !has.long.data)
  }

  # Initialize the log-likelihood values for the various groups (not all may be computed depending
  # on group sizes and arguments)
  g1_longlogl = 0
  g12_survlogl = 0
  g3_logl = 0
  g4_logl = 0

  # Extract the parameter values
  if (args$include.longitudinal) {
    mus=params[args$mu.indexes]
    betas=params[args$beta.indexes]
    if (length(args$psi.indexes)>0) {
      psis=params[args$psi.indexes]
    } else {
      psis = c()
    }
    sigma=params[args$sigma.index]
    nu=params[args$nu.index]
    #eta=params[args$eta.index]
    eta = 0
  }
  if (args$include.survival) {
    alphas = params[args$alpha.indexes]
    if (length(args$phi.indexes)>0) {
      phis=params[args$phi.indexes]
    } else {
      phis = vector()
    }
  }

  # calculate the longitudinal log-likelihood for group 1 (just the longitudinal part)
  if (length(group1) > 0 & args$include.longitudinal) {
    #message("Computing group 1")
    g1_longlogl = group1LongitudinalLogLikelihood(args$long.status[group1,],
        args$num.long.measurements[group1],
        args$long.values[group1,],
        args$long.times[group1,],
        args$long.covariates[group1,,drop=FALSE],
        args$treatment.status[group1],
        long.basis=args$long.basis,
        unique.time.indexes=args$unique.time.indexes, # this is already specific to group 1
        correlation.exponent=args$correlation.exponent,
        mus=mus,betas=betas,psis=psis,sigma=sigma, nu=nu, #eta=eta,
        cache=group1.R.cache)
  }

  # calculate the survival log-likelihood for groups 1 and 2
  if (args$include.survival & (length(group1) + length(group2)) > 0) {
    #message("Computing group 1 and 2")
    sub.surv.basis = NA
    if (!is.na(args$surv.basis)[1]) {
      sub.surv.basis = args$surv.basis[c(group1,group2),,drop=FALSE]
    }
    g12_survlogl = group12SurvivalLogLikelihood(
        covariates=args$surv.covariates[c(group1,group2),,drop=FALSE],
        surv.times=args$surv.times[c(group1,group2)],
        treatment.status=args$treatment.status[c(group1,group2)],
        surv.basis=sub.surv.basis,
        alphas=alphas,
        phis=phis,
        lambdas=lambdas)
  }

  # calculate the log-likelihood for group 3
  if (!args$naive & length(group3) > 0 & args$include.longitudinal & args$include.survival) {
    #message("Computing group 3")
    g3_logl = group3LogLikelihood(long.status=args$long.status[group3,],
        num.long.measurements=args$num.long.measurements[group3],
        long.values=args$long.values[group3,],
        long.times=args$long.times[group3,],
        long.covariates=args$long.covariates[group3,,drop=FALSE],
        surv.covariates=args$surv.covariates[group3,,drop=FALSE],
        treatment.status=args$treatment.status[group3],
        long.knots=args$long.knots,
        surv.basis=args$surv.basis,
        natural.spline=args$natural.spline,
        correlation.exponent=args$correlation.exponent,
        mus=mus,betas=betas,psis=psis,sigma=sigma, nu=nu, #eta=eta,
        cache=group3.R.cache,
        surv.times=args$surv.times,
        surv.times.for.group=args$surv.times[group3],
        censored=args$censored,
        alphas=alphas,
        phis=phis,
        lambdas=lambdas)
  }

  # calculate the log-likelihood for group 4
  if (length(group4) > 0 & args$include.survival) {
    #message("Computing group 4")
    g4_logl = group4LogLikelihood(
        covariates=args$surv.covariates[group4,,drop=FALSE],
        censored=args$censored,
        all.surv.times=args$surv.times,
        surv.times=args$surv.times[group4],
        treatment.status=args$treatment.status[group4],
        surv.basis=args$surv.basis,
        alphas=alphas,
        phis=phis,
        lambdas=lambdas)
  }

  joint_logl = g1_longlogl + g12_survlogl + g3_logl + g4_logl

  # return the negative of the log-likelihood so we can optimize via minimization

  #if (group1.R.cache$eval.num %% 50 == 0) {
  #   message("Joint log-likelihood for eval ", group1.R.cache$eval.num, ": ", joint_logl,
  #       ", longitudinal log-likelihood for group 1: ", g1_longlogl,
  #       ", survival log_likelihood for groups 1 and 2: ", g12_survlogl,
  #       ", log-likelihood for group 3: ", g3_logl,
  #       ", log-likelihood for group 4: ", g4_logl)
  # }
  group1.R.cache$eval.num = group1.R.cache$eval.num + 1

  return (-joint_logl)
}



