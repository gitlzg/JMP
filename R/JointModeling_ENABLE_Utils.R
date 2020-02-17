#
# File: ENABLE_JointModeling_Utils.R
#

library(nlme)
library(mefa) ## for flip the data




#-------------------------------------------------------------------------------------------
# Complete test logic for ENABLE
#-------------------------------------------------------------------------------------------

testENABLE = function(
    use.saved.data = F,
    save.data=F,
    data.file = "JointModeling.RData",
    optimize=F,
    compute.standard.errors=F,
    profile=F,
    naive=F,
    include.longitudinal=T,
    include.survival=T,
    long.knots=NA, # Specific knots, only considered if knots.from.data=F
    surv.knots=NA,# Specific knots, only considered if knots.from.data=F
    knots.from.data=T, # set to true to dynamically compute the specified number of knots from the data
    num.long.knots, # number of long.knots including boundaries if knots.from.data=T; 0 for no spline
    num.surv.knots, # number of surv knots if knots.from.data=T; 0 for no spline
    enabledata.path,
    long.covariate.fields,
    surv.covariate.fields,
    qol.prefix,
    qol.time.prefix,
    has.qol.t0 = F ,
    sample=NA,
    reltol=1e-8,
    outer.iterations=50, # number of outer iterations for calls to constrOptim(); done for all methods except "L-BFGS-B"
    outer.eps=1e-8, # convergence tolerance for constrOptim()
    optim.method="BFGS", #Nelder-Mead",
    enforce.bounds=T, # Controls whether contrOptim() or optim() is called
    natural.spline=F,
    enabledata.table, # set to non-NA value if data file does not need to be read in again
    resample.method=NA, # NA for no resampling or "bootstrap" or "jackknife"
    num.bootstraps=10, # number of bootstrap resampled datasets to generate, only meaningful if resample.method="bootstrap"
    jackknife.size=5, # if resample.method=jackknife, this is the size of the deletion group
    mle=NA,
    init.value.results=NA,
    enable.data_load=NA,
    id.field,
    survival.time.field,
    censoring.status.field,
    treatment.status.field
){

  #-------------------------------------------------------------------------------------------
  # Load ENABLE data
  #-------------------------------------------------------------------------------------------
  enable.data = loadENABLE(
      enabledata.table=enabledata.table,
      enabledata.path = enabledata.path,
      long.covariate.fields = long.covariate.fields,
      surv.covariate.fields = surv.covariate.fields,
      qol.prefix=qol.prefix,
      #num.qol.times=num.qol.times,
      qol.time.prefix=qol.time.prefix,
      has.qol.t0 = has.qol.t0,
      subset.size=NA,
      include.longitudinal=include.longitudinal,
      include.survival=include.survival,
      natural.spline=natural.spline,
      sample=sample,
      id.field = id.field,
      survival.time.field=survival.time.field,
      censoring.status.field=censoring.status.field,
      treatment.status.field=treatment.status.field)
  message("Loaded ENABLE data")
  # print(paste("Number subjects ", enable.data$num.subjects, ", number censored: ", length(which(enable.data$censored)),
  #         ", number of uncensored times: ", length(which(enable.data$uncensored))))
  # print(paste("Number of observed times: ", enable.data$num.events))
  group1 = !enable.data$censored & enable.data$has.qol.data
  group2 = !enable.data$censored & !enable.data$has.qol.data
  group3 = enable.data$censored & enable.data$has.qol.data
  group4 = enable.data$censored & !enable.data$has.qol.data
  # print(paste("Number in group 1 (uncensored, QoL data): ", length(which(group1))))
  # print(paste("Number in group 2 (uncensored, no QoL data): ", length(which(group2))))
  # print(paste("Number in group 3 (censored, QoL data): ", length(which(group3))))
  # print(paste("Number in group 4 (censored, no QoL data): ", length(which(group4))))

  #-------------------------------------------------------------------------------------------
  # Generate resampled datasets via bootstrapping or jackknife if enabled
  #-------------------------------------------------------------------------------------------

  resample = !is.na(resample.method)
  if (resample) {
    bootstrap = (resample.method == "bootstrap")
    if (bootstrap) {
      num.resamples = num.bootstraps
    } else { # jackknife
      num.resamples = floor(enable.data$num.subjects/jackknife.size)
    }
    random.subjects = sample(1:enable.data$num.subjects, enable.data$num.subjects, replace=F)
    resampled.datasets = list()
    jackknife.start = NA
    for (i in 1:num.resamples) {
      message("Computing resampled data set ", i)
      if (bootstrap) {
        resampled.subjects = sample(1:enable.data$num.subjects, enable.data$num.subjects, replace=T)
      } else { # jackknife
        # Determine the range of subjects to delete
        if (is.na(jackknife.start)) {
          jackknife.start = 1
        } else {
          jackknife.start = jackknife.start + jackknife.size
        }
        jackknife.end = jackknife.start + jackknife.size - 1
        # Get indices of all subjects not in deleted group (select from randomly shuffled subject list since
        # treated/untreated may be in blocks)
        resampled.subjects = which(!(random.subjects %in% jackknife.start:jackknife.end))
      }

      resampled.table = enable.data$enabledata.table[resampled.subjects,]
      resampled.datasets[[i]] = loadENABLE(enabledata.table=resampled.table,
          enabledata.path = enabledata.path,
          long.covariate.fields = long.covariate.fields,
          surv.covariate.fields = surv.covariate.fields,
          qol.prefix=qol.prefix,
          #num.qol.times=num.qol.times,
          qol.time.prefix=qol.time.prefix,
          has.qol.t0 = has.qol.t0,
          subset.size=NA,
          include.longitudinal=include.longitudinal,
          include.survival=include.survival,
          natural.spline=natural.spline,
          sample=sample)
    }
  }


  #-------------------------------------------------------------------------------------------
  # Either compute the likelihood once at the initial values or optimize
  #-------------------------------------------------------------------------------------------

  if ( use.saved.data) {
    message("Loading saved data...")
    #load(data.file)
    # Set the param indexes from the loaded saved data
    param.indexes = getParamIndexes(
      num.long.fixed.effects=ncol(enable.data_load$long.covariates),
      num.surv.fixed.effects=ncol(enable.data_load$surv.covariates),
      num.long.knots=enable.data_load$num.long.knots,
      num.surv.knots=enable.data_load$num.surv.knots,
      num.events=enable.data_load$num.events,
      include.survival=enable.data_load$include.survival,
      include.longitudinal=enable.data_load$include.longitudinal,
      natural.spline=enable.data_load$natural.spline)
    long.knots=enable.data_load$long.knots
    surv.knots=enable.data_load$surv.knots
    message("...finished loading saved data ")
  } else {
    message("Did not find saved R data ", data.file, ", computing from scratch")

    #-------------------------------------------------------------------------------------------
    # Process spline knot parameters
    #-------------------------------------------------------------------------------------------

    if (knots.from.data) { # Dynamically compute the knots
      knots = computeSplineKnots(enable.data=enable.data, num.long.knots=num.long.knots, num.surv.knots=num.surv.knots)
      long.knots = knots$long.knots
      surv.knots = knots$surv.knots
    } else { # Use specific knots
      num.long.knots = length(long.knots)
      num.surv.knots = length(surv.knots)
    }

    # add the knot settings to the enable.data structure
    enable.data$long.knots=long.knots
    enable.data$num.long.knots=num.long.knots
    enable.data$surv.knots=surv.knots
    enable.data$num.surv.knots = num.surv.knots

    #-------------------------------------------------------------------------------------------
    # Get likelihood parameter indexes
    #-------------------------------------------------------------------------------------------

    param.indexes = getParamIndexes(
        num.long.fixed.effects=ncol(enable.data$long.covariates),
        num.surv.fixed.effects=ncol(enable.data$surv.covariates),
        num.long.knots=num.long.knots,
        num.surv.knots=num.surv.knots,
        num.events=enable.data$num.events,
        include.survival=enable.data$include.survival,
        include.longitudinal=enable.data$include.longitudinal,
        natural.spline=natural.spline)

    # create the parameter initial values
    init.value.results = getENABLEInitialValues(param.indexes, enable.data=enable.data,
                                                treatment.status.field=treatment.status.field,
                                                id.field = id.field,
                                                survival.time.field=survival.time.field,
                                                censoring.status.field=censoring.status.field)
    param.initial.values = init.value.results$init.param.values
    message("Initial values: ", paste(param.initial.values, collapse=", "))
    if (optimize) {
      optim.method=optim.method
      control = list(#trace=6,# give some tracing
          #ndeps=1e-2,
          #maxit=10000, # control the maximum number of iterations
          reltol=reltol # relative tolerance for convergence
      #alpha=1,# Nelder-Mead reflection factor (default is 1)
      #beta=.75, # Nelder-Mead contraction factor (default is .5); increase to not let the polytope get too small
      #gamma=3, # Nelder-Mead expansion factor (default is 2); increase to go farther in good directions
      #REPORT=10
      ) # tracing frequency in iterations
      # set the scaling for the parameters to the absolute value of the initial value
      control$parscale=abs(param.initial.values)
      # create the parameter bounds
      param.bounds = getENABLEParameterBounds(param.indexes, param.initial.values)
      message("Lower bounds: ", paste(param.bounds$param.lower.bounds, collapse=", "))
      message("Upper bounds: ", paste(param.bounds$param.upper.bounds, collapse=", "))
    } else {
      optim.method=NA
      param.bounds = list(param.lower.bounds=NA, param.upper.bounds=NA)
      control = NA
    }
    start = proc.time()[3]

    # compute or optimize the log-likelihood
    mle = optimizeJointLikelihood(long.values = enable.data$long.values,
        long.times = enable.data$long.times,
        surv.times = enable.data$surv.times,
        censored=enable.data$censored,
        treatment.status=enable.data$treatment.status,
        long.covariates=enable.data$long.covariates,
        surv.covariates=enable.data$surv.covariates,
        long.knots=long.knots,
        surv.knots=surv.knots,
        natural.spline=natural.spline,
        correlation.exponent=2,
        param.init.values=param.initial.values,
        param.lower.bounds=param.bounds$param.lower.bounds,
        param.upper.bounds=param.bounds$param.upper.bounds,
        optim.method=optim.method,
        control=control,
        naive=naive,
        outer.iterations=outer.iterations,
        outer.eps=outer.eps,
        enforce.bounds=enforce.bounds,
        include.survival=enable.data$include.survival,
        include.longitudinal=enable.data$include.longitudinal,
        compute.standard.errors=compute.standard.errors)

    # add the knots
    mle$long.knots = long.knots
    mle$surv.knots = surv.knots

    total.time = proc.time()[3] - start
    #print(paste("Total elapsed time: ", total.time, " sec"))
    if (profile & optimize) {
      print(summaryRprof("JointModelingProf.out"))
    }

    if (resample) {
      resampled.mles = list()
      success.index = 0
      for (i in 1:num.resamples) {

        message("Performing resample test ", i)
        resampled.data = resampled.datasets[[i]]
        if (knots.from.data) { # Dynamically compute the knots
          knots = computeSplineKnots(enable.data=resampled.data, num.long.knots=num.long.knots, num.surv.knots=num.surv.knots)
          long.knots = knots$long.knots
          surv.knots = knots$surv.knots
        } else { # Use specific knots
          num.long.knots = length(long.knots)
          num.surv.knots = length(surv.knots)
        }

        # add the knot settings to the enable.data structure
        resampled.data$long.knots=long.knots
        resampled.data$num.long.knots=num.long.knots
        resampled.data$surv.knots=surv.knots
        resampled.data$num.surv.knots = num.surv.knots

        # create the parameter initial values for the resampled data set
        resampled.init.value.results = try(getENABLEInitialValues(param.indexes, enable.data=resampled.data,
                                                                  treatment.status.field=treatment.status.field,
                                                                  id.field = id.field,
                                                                  survival.time.field=survival.time.field,
                                                                  censoring.status.field=censoring.status.field))
        if (inherits(resampled.init.value.results, "try-error")) {
          message("Error computing initial values for resampled data ", i)
          # In case needed, change warnings back to normal
          options(warn=1)
          next
        }

        resampled.param.initial.values = resampled.init.value.results$init.param.values
        resampled.param.bounds = getENABLEParameterBounds(param.indexes, resampled.param.initial.values)
        control$parscale=abs(resampled.param.initial.values)

        # compute or optimize the log-likelihood
        resampled.mle = try(optimizeJointLikelihood(long.values = resampled.data$long.values,
                long.times = resampled.data$long.times,
                surv.times = resampled.data$surv.times,
                censored=resampled.data$censored,
                treatment.status=resampled.data$treatment.status,
                long.covariates=resampled.data$long.covariates,
                surv.covariates=resampled.data$surv.covariates,
                long.knots=long.knots,
                surv.knots=surv.knots,
                natural.spline=natural.spline,
                correlation.exponent=2,
                param.init.values=resampled.param.initial.values,
                param.lower.bounds=resampled.param.bounds$param.lower.bounds,
                param.upper.bounds=resampled.param.bounds$param.upper.bounds,
                optim.method=optim.method,
                control=control,
                naive=naive,
                outer.iterations=outer.iterations,
                outer.eps=outer.eps,
                enforce.bounds=enforce.bounds,
                include.survival=resampled.data$include.survival,
                include.longitudinal=resampled.data$include.longitudinal,
                compute.standard.errors=F))
        if (!inherits(resampled.mle, "try-error")) {
          success.index = success.index + 1
          resampled.mles[[success.index]] = resampled.mle
        } else {
          message("Error optimizing for resampled data ", i)
          # In case needed, change warnings back to normal
          options(warn=1)
        }
      }
      resampled.estimates = matrix(0, nrow=success.index, ncol=length(param.indexes$param.names))
      for (i in 1:success.index) {
        resampled.estimates[i,] = resampled.mles[[i]]$optim.results$par
      }

      # save the resampled estimates
      mle$resampled.estimates=resampled.estimates

      # Use resampled values to approximate the estimate covariance matrix
      if (bootstrap) {
        mle$resampled.cov.mat = cov(resampled.estimates)
      } else {
        k = nrow(resampled.estimates)
        mle$resampled.cov.mat =((k-1)^2/k)*cov(resampled.estimates)
      }
      mle$resampled.standard.errors = sqrt(diag(mle$resampled.cov.mat))

      # Compute resample bias-corrected estimate
      resampled.estimate.means = apply(resampled.estimates, 2, mean)
      if (bootstrap) {
        mle$resampled.par = 2*mle$optim.results$par - resampled.estimate.means
      } else {
        k = nrow(resampled.estimates)
        mle$resampled.par = mle$optim.results$par - (k-1)*(resampled.estimate.means-mle$optim.results$par)
      }
    }

    if (save.data) {
      save(mle, init.value.results, enable.data, total.time, file=data.file)
    }
    enable.data_load=enable.data
  }

  # print(mle)
  # print(init.value.results)

  results = list(enable.data=enable.data_load, mle=mle,
      init.value.results=init.value.results,
      long.knots=long.knots, surv.knots=surv.knots, param.indexes=param.indexes)
  return (results)
}

#-------------------------------------------------------------------
# BIC-based selection for the number of survival and longitudinal knots
#-------------------------------------------------------------------

knotSelectionViaBIC = function(
    optimize=F,
    naive=F,
    include.longitudinal=T,
    include.survival=T,
    long.knot.nums,
    surv.knot.nums,
    enabledata.path,
    enabledata,
    long.covariate.fields,
    surv.covariate.fields,
    qol.prefix,
    qol.time.prefix,
    has.qol.t0 = F ,
    sample=NA,
    reltol=1e-8,
    outer.iterations=50, # number of outer iterations for calls to constrOptim(); done for all methods except "L-BFGS-B"
    outer.eps=1e-8, # convergence tolerance for constrOptim()
    optim.method="BFGS", #Nelder-Mead",
    enforce.bounds=T, # Controls whether contrOptim() or optim() is called
    natural.spline=F,
    id.field,
    survival.time.field,
    censoring.status.field,
    treatment.status.field
){
  knot.combinations = expand.grid(num.long.knots=long.knot.nums, num.surv.knots=surv.knot.nums)
  knot.combinations$BIC = rep(NA, nrow(knot.combinations))
  knot.combinations$AIC = rep(NA, nrow(knot.combinations))
  enabledata.table=enabledata
  for (i in 1:nrow(knot.combinations)) {
    num.long.knots = knot.combinations$num.long.knots[i]
    num.surv.knots = knot.combinations$num.surv.knots[i]
    message("Testing AIC/BIC for ", num.long.knots, " long knots and ", num.surv.knots, " surv knots...")
    result = try(testENABLE(
            use.saved.data = F,
            save.data=F,
            optimize=optimize,
            compute.standard.errors=F,
            enabledata.table=enabledata.table,
            profile=F,
            naive=naive,
            include.longitudinal=include.longitudinal,
            include.survival=include.survival,
            long.knots=NA, surv.knots=NA,
            knots.from.data=T,
            num.long.knots=num.long.knots,
            num.surv.knots=num.surv.knots,
            enabledata.path=enabledata.path,
            long.covariate.fields =long.covariate.fields,
            surv.covariate.fields =surv.covariate.fields,
            qol.prefix=qol.prefix,
            #num.qol.times=num.qol.times,
            qol.time.prefix=qol.time.prefix,
            has.qol.t0=has.qol.t0,
            sample=sample,
            reltol=reltol,
            outer.iterations=outer.iterations,
            outer.eps=outer.eps,
            optim.method=optim.method,
            enforce.bounds=enforce.bounds,
            natural.spline=natural.spline,
            id.field = id.field,
            survival.time.field=survival.time.field,
            censoring.status.field=censoring.status.field,
            treatment.status.field=treatment.status.field
            ))
    if (!inherits(result, "try-error")) {
      enabledata.table=result$enable.data$enabledata.table
      num.params = length(result$param.indexes$param.names)
      num.samples = result$enable.data$num.subjects
      neg.log.likelihood = result$mle$neg.log.likelihood
      bic = 2*neg.log.likelihood + log(num.samples)*num.params
      aic = 2*neg.log.likelihood + 2*num.params
      message("AIC/BIC for ", num.long.knots, " long knots and ", num.surv.knots, " surv knots, ", num.params, " params, and ", num.samples, " samples: ", aic, "/", bic)
      if (!is.infinite(aic) & !is.infinite(bic)) {
        knot.combinations$BIC[i] = bic
        knot.combinations$AIC[i] = aic
      } else {
        message("Infinite AIC/BIC, not using")
      }
    } else {
      message("Error computing AIC/BIC for ", num.long.knots, " long knots and ", num.surv.knots, " surv knots")
      enabledata.table=NA
      # Change warnings back to normal
      options(warn=1)
    }
  }
  #print("knot.combinations:")
  #print(knot.combinations)
  return (knot.combinations)
}


#-------------------------------------------------------------------
# Compute spline knots
#-------------------------------------------------------------------
computeSplineKnots = function(enable.data, num.surv.knots, num.long.knots) {

  message("Computing spline knots")

  if (num.surv.knots > 0) {
    probs = seq(from=0,to=1,length.out=num.surv.knots)
    surv.knots = quantile(enable.data$surv.times, probs)
    message("surv.knots: ", paste(surv.knots, collapse=","))
  } else {
    surv.knots = c()
  }

  if (num.long.knots > 0) {
    probs = seq(from=0,to=1,length.out=num.long.knots)
    long.times = as.vector(enable.data$long.times)
    pos.long.times = long.times[which(long.times >= 0)]
    long.knots = quantile(pos.long.times, probs)
    message("long.knots: ", paste(long.knots, collapse=","))
  } else {
    long.knots = c()
  }

  results = list()
  results$long.knots=long.knots
  results$surv.knots=surv.knots
  return (results)
}

#-------------------------------------------------------------------
# Loads the ENABLE data
#-------------------------------------------------------------------
loadENABLE = function(
    enabledata.table,
    subset.size = NA, #100 # NA to process all of the data, number of process only the first X subjects
    long.covariate.fields, # name of column headers hold covariates; leave empty for no covariates
    surv.covariate.fields, # name of column headers hold covariates; leave empty for no covariates
    include.survival=T,
    include.longitudinal=T,
    natural.spline=F,
    enabledata.path, #"../data/enable2_wide_rob.csv"
    qol.prefix, #"FACTPAL_"
    qol.time.prefix, #"t_"
    has.qol.t0 = F, #T
    sample=NA,
    id.field,
    survival.time.field,
    censoring.status.field,
    treatment.status.field) {

  #-------------------------------------------------------------------
  # Load the data
  #-------------------------------------------------------------------

  if (sum(!is.na(enabledata.table))==0) {
    message("Reading ", enabledata.path, "...")
    enabledata.table = read.table(enabledata.path, header=T, sep=",")
    # remove baseline missing values
    enabledata.table=enabledata.table[!(NA%in%enabledata.table[,long.covariate.fields]),]
    enabledata.table=enabledata.table[!(NA%in%enabledata.table[,surv.covariate.fields]),]
    message("Subjects missing baseline covariates values (if any) are removed in the analysis.")

    # whether simulation
    if (!is.na(sample)) {
      enabledata.table = subset(enabledata.table, rep == sample)
      enabledata.table = enabledata.table[,2:ncol(enabledata.table)]
    }
    message("...finished reading ", enabledata.path)
  }
  message("Number of subjects: ", nrow(enabledata.table))
  message("Number of fields: ", ncol(enabledata.table))
  message("Names of ENABLE data columns: ")
  message(paste(names(enabledata.table), collapse=","))
  if (!is.na(subset.size)) {
    message("Subsetting to first ", subset.size, " subjects...")
    enabledata.table = enabledata.table[1:subset.size,]
  }
  enabledata.table = enabledata.table

  #-------------------------------------------------------------------
  # Extract and compute some useful values from the input data
  #-------------------------------------------------------------------

  treatment.status = enabledata.table[,treatment.status.field]

  long.covariates = as.matrix(apply(enabledata.table[,long.covariate.fields, drop=F], c(1,2), as.numeric))
  surv.covariates = as.matrix(apply(enabledata.table[,surv.covariate.fields, drop=F], c(1,2), as.numeric))
  surv.times = enabledata.table[,survival.time.field]
  censored = enabledata.table[,censoring.status.field] == 0
  num.censored = length(which(censored))
  uncensored = enabledata.table[,censoring.status.field] == 1
  num.events = length(unique(surv.times))
  num.subjects = nrow(enabledata.table)

  #-------------------------------------------------------------------
  # Compute:
  # - Overall proportion of censored subjects
  # - Proportion of censored subjects after the last real death time
  #-------------------------------------------------------------------

  message("Number of subjects ", num.subjects, ", num censored ", num.censored)
  overall.censored.proportion = num.censored/num.subjects
  last.death.time = max(surv.times[uncensored])
  num.censored.after.last.death = length(which(surv.times > last.death.time))
  proportion.censored.after.last.death = num.censored.after.last.death/num.censored
  message("Proportion of censored among all subjects: ", overall.censored.proportion)
  message("Number censored after last death time: ", num.censored.after.last.death)
  message("Proportion of censored after last death time: ", proportion.censored.after.last.death)

  #-------------------------------------------------------------------
  # IMPORTANT: make sure the last event time is a death (so, if censored, change to death);
  # this ensures that all censoring events happen before a death
  #-------------------------------------------------------------------
  surv.times = enabledata.table[,survival.time.field]
  last.event.time = max(surv.times)
  last.event = tail(which(surv.times == last.event.time), n=1)
  if (enabledata.table[last.event, censoring.status.field] == 0) {
    #stop("Last event time was censored, halting processing")
    message("Last event was censored, changing to death")
    enabledata.table[last.event, censoring.status.field] = 1
    # add 0.1 to the last observed time to avoid death time ties
    enabledata.table[last.event, survival.time.field] = enabledata.table[last.event, survival.time.field]+0.1
    surv.times[last.event]=enabledata.table[last.event, survival.time.field]
    # recompute the censored and uncensored counts
    censored = enabledata.table[,censoring.status.field] == 0
    uncensored = enabledata.table[,censoring.status.field] == 1
  } else {
    message("Last event was a death, no need to change")
  }
  #print("length(unique(surv.times[uncensored])):")
  #print(length(unique(surv.times[uncensored])))

  if (include.longitudinal) {
    qol.cols = sapply(colnames(enabledata.table), function(x) {grep(qol.prefix, x)})
    num.qol.times = length(which(qol.cols == 1))-1 # need to take off qol time 0

    qol.col.names = paste(qol.prefix, 1:num.qol.times, sep="")
    if (has.qol.t0) {
      qol.time.col.names = paste(qol.time.prefix, 0:(num.qol.times-1), sep="")
    } else {
      qol.time.col.names = paste(qol.time.prefix, 1:num.qol.times, sep="")
    }

    minDiff=min(surv.times-enabledata.table[,qol.time.col.names],na.rm=T)
    if(minDiff<=0)stop("Some follow up time is larger than Survival time, double check the calcuated times.")
    minSurv=min(surv.times)
    if(minSurv<=0)stop("Some survival time is negative, double check the calcuated times.")
    minTimes=min(enabledata.table[,qol.time.col.names],na.rm=T)
    if(minTimes<=0)stop("Some follow up time is negative, double check the calcuated times.")

    #message("QoL col names: ", paste(qol.col.names, collapse=", "))
    #message("Old QoL col names: ", paste(colnames(enabledata.table)[4:(3+num.qol.times)], collapse=", "))
    #message("QoL time col names: ", paste(qol.time.col.names, collapse=", "))

    qol.values = enabledata.table[,qol.col.names]
    #qol.values = enabledata.table[,4:(3+num.qol.times)]
    #qol.time.start = 9+num.qol.times
    #qol.time.end = 8+2*num.qol.times
    #if (!has.qol.t0) {
    #  qol.time.start = qol.time.start-1
    #  qol.time.end = qol.time.end-1
    #}
    #qol.times = enabledata.table[,qol.time.start:qol.time.end]

    qol.times = enabledata.table[,qol.time.col.names]
    # turn the times in to time before death
    qol.times = surv.times-qol.times
    qol.status = apply(qol.values, c(1,2), function(x) {as.numeric(!is.na(x))})
    num.qol.measurements = apply(qol.status, 1, sum)
  } else {
    qol.values = NA
    qol.times = NA
    qol.status = NA
    num.qol.measurements=0
  }

  results = list()
  results$enabledata.table = enabledata.table
  results$num.subjects = num.subjects
  results$treatment.status = as.vector(treatment.status)
  results$long.covariates = long.covariates
  results$surv.covariates = surv.covariates
  results$surv.times = as.matrix(surv.times)
  results$long.values = as.matrix(qol.values)
  results$long.times = as.matrix(qol.times)
  results$censored = as.vector(censored)
  results$uncensored = as.vector(uncensored)
  results$num.censored.after.last.death = num.censored.after.last.death
  results$has.qol.data = num.qol.measurements > 0
  results$treatment.status = as.vector(treatment.status)
  results$num.events = num.events
  results$include.longitudinal = include.longitudinal
  results$include.survival = include.survival
  results$long.covariate.fields = long.covariate.fields
  results$surv.covariate.fields = surv.covariate.fields
  results$natural.spline=natural.spline

  return (results)
}

#-------------------------------------------------------------------
# Compute initial values for all of the parameters
#-------------------------------------------------------------------
getENABLEInitialValues = function(param.indexes, enable.data,
                                  treatment.status.field,
                                  id.field,
                                  survival.time.field,
                                  censoring.status.field) {

  init.param.values = rep(0, length(param.indexes$param.names))
  names(init.param.values) = param.indexes$param.names
  init.param.ses = rep(0, length(param.indexes$param.names))
  names(init.param.ses) = param.indexes$param.names
  mu.cov.mat = NA
  beta.cov.mat = NA
  alpha.cov.mat = NA

  if (param.indexes$include.longitudinal) {
    init.values = computeLongitudinalInitialValues(enable.data,
                                                   treatment.status.field=treatment.status.field,
                                                   id.field = id.field,
                                                   survival.time.field=survival.time.field,
                                                   censoring.status.field=censoring.status.field)
    mu.cov.mat = init.values$mu.cov.mat
    beta.cov.mat = init.values$beta.cov.mat
    init.param.values[param.indexes$mu.indexes] = init.values$mu
    init.param.values[param.indexes$beta.indexes] = init.values$beta
    init.param.values[param.indexes$sigma.index] = init.values$sigma
    init.param.values[param.indexes$nu.index] = init.values$nu
    #init.param.values[param.indexes$eta.index] = init.values$eta
    init.param.values[param.indexes$psi.indexes]= init.values$psi

    init.param.ses[param.indexes$mu.indexes] = init.values$mu.se
    init.param.ses[param.indexes$beta.indexes] = init.values$beta.se
    init.param.ses[param.indexes$sigma.index] = init.values$sigma.se
    init.param.ses[param.indexes$nu.index] = init.values$nu.se
    #init.param.values[param.indexes$eta.index] = init.values$eta
    init.param.ses[param.indexes$psi.indexes]= init.values$psi.se

  }

  if (param.indexes$include.survival) {
    init.values = computeSurvivalInitialValues(enable.data,
                                               id.field = id.field,
                                               survival.time.field=survival.time.field,
                                               censoring.status.field=censoring.status.field,
                                               treatment.status.field=treatment.status.field)
    alpha.cov.mat = init.values$alpha.cov.mat
    init.param.values[param.indexes$alpha.indexes] = init.values$alpha
    init.param.values[param.indexes$phi.indexes] = init.values$phi

    init.param.ses[param.indexes$alpha.indexes] = init.values$alpha.se
    init.param.ses[param.indexes$phi.indexes] = init.values$phi.se
  }

  result = list()
  result$init.param.values = init.param.values
  result$init.param.ses = init.param.ses
  result$mu.cov.mat = mu.cov.mat
  result$beta.cov.mat = beta.cov.mat
  result$alpha.cov.mat = alpha.cov.mat

  return (result)
}

#-------------------------------------------------------------------
# Compute initial values for longitudinal parameters based on lme()
#-------------------------------------------------------------------
computeLongitudinalInitialValues = function(enable.data,
                                            group13.indep=F,
                                            treatment.status.field,
                                            id.field,
                                            survival.time.field,
                                            censoring.status.field) {

  message("Computing longitudinal initial values using LME")

  # compute the longitudinal model initial values

  if (enable.data$num.long.knots == 0 & group13.indep) {
    # Combine groups 1 and 3 if no spline
    group1 = enable.data$has.qol.data
    g1_e3=enable.data$enabledata.table[which(group1),]
  } else {
    # Use just group 1
    group1 = !enable.data$censored & enable.data$has.qol.data
    g1_e3=enable.data$enabledata.table[which(group1),]
  }

  # Use just group 3
  #group1 = enable.data$censored & enable.data$has.qol.data
  #g1_e3=enable.data$enabledata.table[which(group1),]

  # check if all of the subjects in group 1 have the same gender
#  gp_gender = g1_e3$sex
#  num.females = length(which(gp_gender == 0))
#  include.sex = T
#  if (num.females == 0 | num.females == nrow(g1_e3)) {
#    # all female or all male
#    message("All group 1 members are female or male, not estimating sex for lme")
#		include.sex = F
#  }

  ############################## flip group 1 data

  long_gp1=rep(g1_e3, each=ncol(enable.data$long.times))
  tstar_vars=colnames(enable.data$long.times)
  qol_vars=colnames(enable.data$long.values)
  tstar_m=g1_e3[tstar_vars]
  qol_m=g1_e3[qol_vars]

  long_gp1$tstar=c(t(tstar_m))
  long_gp1$qol=c(t(qol_m))

  ############################### create the retrospective time scale
  long_gp1$t=long_gp1$tstar

  long_gp1$t=(long_gp1$t>0)*(long_gp1$t)

  ############################### create B-spline basis functions for the control and treatment groups
  B_base = getSplineBasis(knots=enable.data$long.knots, times=long_gp1$t, natural.spline=enable.data$natural.spline)

  if (enable.data$num.long.knots > 0) {
    if (enable.data$natural.spline) {
      num.basis.params = enable.data$num.long.knots
    } else {
      num.basis.params = enable.data$num.long.knots + 2
    }
  } else {
    num.basis.params = 1
  }

  formula.str = "qol~"
  for (i in 1:num.basis.params) {
    mu.name = paste("mu_b", i, sep="")
    long_gp1[[mu.name]] = B_base[,i]
    formula.str = paste(formula.str, mu.name, "+", sep="")
  }
  for (i in 1:num.basis.params) {
    beta.name = paste("beta_b", i, sep="")
    long_gp1[[beta.name]] = long_gp1[,treatment.status.field]*B_base[,i]
    formula.str = paste(formula.str, beta.name, "+", sep="")
  }
  if (length(enable.data$long.covariate.fields) > 0) {
    for (i in 1:length(enable.data$long.covariate.fields)) {
      cov.name = enable.data$long.covariate.fields[i]
      #if (cov.name != "sex" | include.sex) {
      formula.str = paste(formula.str, cov.name, sep="")
      if (i < length(enable.data$long.covariate.fields)) {
        formula.str = paste(formula.str, "+", sep="")
      }
    }
  }
  formula.str = paste(formula.str, "-1", sep="")

  message("LME formula for longitudinal initial values: ", formula.str)

  random.formula.str=as.formula(paste("~1|",id.field,sep=""))

  # make warnings appears as errors for duration of call
  options(warn=2)

  g1mix=lme(as.formula(formula.str),
#        random=~1|id,data=long_gp1,correlation=corGaus(form=~t|id),method="ML",na.action="na.omit")
      random=random.formula.str,data=long_gp1,method="ML",na.action="na.omit")

  # Change warnings back to normal
  options(warn=1)

  ################# extract the parameter estimates

  beta=g1mix$coefficients$fixed
  beta.SE = sqrt(diag(g1mix$varFix))
  message("Fixed coef estimates: ", paste(beta, collapse=", "))
  message("Fixed coef SEs: ", paste(beta.SE, collapse=", "))
  resid_sd=VarCorr(g1mix)[,2]
  eta=1/coef(g1mix$modelStruct$corStruct,unconstrained=FALSE)^2
  #eta=g1mix$modelStruct$corStruct^2
  para=c(beta,resid_sd,eta)
  para.SE=c(beta.SE,0,0)
  para=sapply(para, as.numeric)
  para.SE=sapply(para.SE, as.numeric)

  # save initial values is result list
  results = list()
  results$mu = para[1:num.basis.params]
  results$mu.se = para.SE[1:num.basis.params]
  results$mu.cov.mat = g1mix$varFix[1:num.basis.params, 1:num.basis.params]
  beta.range = (num.basis.params+1):(2*num.basis.params)
  results$beta = para[beta.range]
  results$beta.se = para.SE[beta.range]
  results$beta.cov.mat = g1mix$varFix[beta.range, beta.range]
  if (length(enable.data$long.covariate.fields) > 0) {
    for (i in 1:length(enable.data$long.covariate.fields)) {
      cov.name = enable.data$long.covariate.fields[i]
      message("Getting initial value of ", cov.name)
      results$psi[i] = para[(2*num.basis.params+i)]
      results$psi.se[i] = para.SE[(2*num.basis.params+i)]
    }
    #results$psi = para[(2*num.basis.params+1):(2*num.basis.params+length(enable.data$covariate.fields))]
    offset = length(enable.data$long.covariate.fields) +1
  } else {
    offset = 1
  }

  # HACK: Hard code psi2 and psi4
  #results$psi[2] = 5
  #results$psi[4] = 4

  message("Initial psi values: ", paste(results$psi, collapse=", "))

  results$sigma = para[(2*num.basis.params+offset)]
  results$sigma.se = para.SE[(2*num.basis.params+offset)]
  results$nu = para[(2*num.basis.params+offset+1)]
  results$nu.se = para.SE[(2*num.basis.params+offset+1)]
  #results$eta = para[(2*num.basis.params+offset+2)]

  message("Finished computing longitudinal initial values.")

  return (results)
}

#-------------------------------------------------------------------
# Fit a Cox proportional hazards model to get initial values for the survival model
#-------------------------------------------------------------------
computeSurvivalInitialValues = function(enable.data,
                                        id.field,
                                        survival.time.field,
                                        censoring.status.field,
                                        treatment.status.field) {

  init.values = list()

  message("Computing time-varying covariates Cox model")

  ############################### transform the data into a time-dependent format
  uniq.event.times=unique(sort(enable.data$enabledata.table[,survival.time.field][enable.data$enabledata.table[,censoring.status.field]==1]))

  n=nrow(enable.data$enabledata.table)
  message("n: ", n)
  reps=rep(NA,n)

  for (i in 1:n){
    reps[i]=sum(as.numeric(enable.data$enabledata.table[,survival.time.field][i]>uniq.event.times))+1
  }

  long.data=rep(enable.data$enabledata.table, reps)
  cumu=c(0,cumsum(reps))

  long.data$tstart=rep(NA,cumu[n+1])
  long.data$tstop=rep(NA,cumu[n+1])
  long.data$long.censor=rep(0,cumu[n+1])

  for (i in 1:n){
    if(reps[i]>1){
      long.data$tstart[(cumu[i]+1):cumu[i+1]]=c(0,uniq.event.times[1:(reps[i]-1)])
      }else{
        long.data$tstart[(cumu[i]+1):cumu[i+1]]=c(0)
        }
    long.data$tstop[(cumu[i]+1):cumu[i+1]]=uniq.event.times[1:reps[i]]
    long.data$tstop[cumu[i+1]]=enable.data$enabledata.table[,survival.time.field][i]
    long.data$long.censor[cumu[i+1]]=enable.data$enabledata.table[,censoring.status.field][i]
  }

  B_base = getSplineBasis(knots=enable.data$surv.knots, times=long.data$tstop, natural.spline=enable.data$natural.spline)

  if (enable.data$num.surv.knots > 0) {
    if (enable.data$natural.spline) {
      num.basis.params = enable.data$num.surv.knots
    } else {
      num.basis.params = enable.data$num.surv.knots +2
    }
  } else {
    num.basis.params = 1
  }

  formula.str = "Surv(tstart, tstop, long.censor) ~ "
  for (i in 1:num.basis.params) {
    alpha.name = paste("alpha_b", i, sep="")
    long.data[[alpha.name]] = long.data[,treatment.status.field]*B_base[,i]
    formula.str = paste(formula.str, alpha.name, sep="")
    if (i < num.basis.params) {
      formula.str = paste(formula.str, "+", sep="")
    }
  }
  if (length(enable.data$surv.covariate.fields) > 0) {
    for (i in 1:length(enable.data$surv.covariate.fields)) {
      formula.str = paste(formula.str, "+", enable.data$surv.covariate.fields[i], sep="")
    }
  }

  message("Cox formula for survival initial values: ", formula.str)

  # make warnings appears as errors for duration of call
  options(warn=2)

  coxph.results=coxph(formula = as.formula(formula.str), data=long.data)
  coxph.SEs = sqrt(diag(coxph.results$var))
  message("Coxph SEs: ", paste(coxph.SEs, collapse=", "))

  # Change warnings back to normal
  options(warn=1)

  init.values$alpha.cov.mat = coxph.results$var[1:num.basis.params, 1:num.basis.params]
  init.values$alpha = coxph.results$coefficients[1:num.basis.params]
  init.values$alpha.se = coxph.SEs[1:num.basis.params]
  init.values$phi = coxph.results$coefficients[(num.basis.params+1):(num.basis.params+length(enable.data$surv.covariate.fields))]
  init.values$phi.se = coxph.SEs[(num.basis.params+1):(num.basis.params+length(enable.data$surv.covariate.fields))]

  return (init.values)
}

#-------------------------------------------------------------------
# Get bounds for all of the parameters, edit as desired
#-------------------------------------------------------------------

getENABLEParameterBounds = function(param.indexes, param.init.values, prop=4) {
  return (getInitialValueBounds(param.indexes, param.init.values, prop=prop))
  #return (getOldBounds(param.indexes))
}

getLowerBound = function(init.value, proportion, min.bound=NA) {
  lower.bound = init.value - abs(init.value*proportion)
  if (!is.na(min.bound)) {
    lower.bound = sapply(lower.bound, function(x) {
          if (x < min.bound) {
            return (min.bound)
          } else {
            return (x)
          }
        })
  }
  return (lower.bound)
}
getUpperBound = function(init.value, proportion) {
  upper.bound = init.value + abs(init.value*proportion)
  return (upper.bound)
}

getInitialValueBounds = function(param.indexes, param.init.values, prop) {

  param.lower.bounds = rep(0, length(param.indexes$param.names))
  param.upper.bounds = rep(0, length(param.indexes$param.names))
  names(param.lower.bounds) = param.indexes$param.names
  names(param.upper.bounds) = param.indexes$param.names

  # set initial values and bounds for longitudinal model parameters
  if (param.indexes$include.longitudinal) {
    param.lower.bounds[param.indexes$mu.indexes] = getLowerBound(param.init.values[param.indexes$mu.indexes], prop)
    param.upper.bounds[param.indexes$mu.indexes] = getUpperBound(param.init.values[param.indexes$mu.indexes], prop)
    param.lower.bounds[param.indexes$beta.indexes] = getLowerBound(param.init.values[param.indexes$beta.indexes], prop)
    param.upper.bounds[param.indexes$beta.indexes] = getUpperBound(param.init.values[param.indexes$beta.indexes], prop)
    if (length(param.indexes$psi.indexes) > 0) {
      # order is age, sex, FACTPal_0
      param.lower.bounds[param.indexes$psi.indexes] = getLowerBound(param.init.values[param.indexes$psi.indexes],prop)
      param.upper.bounds[param.indexes$psi.indexes] = getUpperBound(param.init.values[param.indexes$psi.indexes],prop)
    }
    param.lower.bounds[param.indexes$sigma.index] = getLowerBound(param.init.values[param.indexes$sigma.index],prop, min.bound=0.00000001)
    param.upper.bounds[param.indexes$sigma.index] = getUpperBound(param.init.values[param.indexes$sigma.index],prop)
    param.lower.bounds[param.indexes$nu.index] = getLowerBound(param.init.values[param.indexes$nu.index],prop, min.bound=0.00000001)
    param.upper.bounds[param.indexes$nu.index] = getUpperBound(param.init.values[param.indexes$nu.index],prop)
    #param.lower.bounds[param.indexes$eta.index] = 0.00000001
    #param.upper.bounds[param.indexes$eta.index] = 25
  }

  # set initial values and bounds for survival model parameters
  if (param.indexes$include.survival) {
    param.lower.bounds[param.indexes$alpha.indexes] = getLowerBound(param.init.values[param.indexes$alpha.indexes],prop)
    param.upper.bounds[param.indexes$alpha.indexes] = getUpperBound(param.init.values[param.indexes$alpha.indexes],prop)
    if (length(param.indexes$phi.indexes) > 0) {
      param.lower.bounds[param.indexes$phi.indexes] = getLowerBound(param.init.values[param.indexes$phi.indexes],prop)
      param.upper.bounds[param.indexes$phi.indexes] = getUpperBound(param.init.values[param.indexes$phi.indexes],prop)
    }
  }

  results = list()
  results$param.lower.bounds = param.lower.bounds
  results$param.upper.bounds = param.upper.bounds

  return (results)
}

getOldBounds = function(param.indexes) {

  param.lower.bounds = rep(0, length(param.indexes$param.names))
  param.upper.bounds = rep(0, length(param.indexes$param.names))
  names(param.lower.bounds) = param.indexes$param.names
  names(param.upper.bounds) = param.indexes$param.names

  # set initial values and bounds for longitudinal model parameters
  if (param.indexes$include.longitudinal) {
    param.lower.bounds[param.indexes$mu.indexes] = -10^9
    param.upper.bounds[param.indexes$mu.indexes] = 10^9
    param.lower.bounds[param.indexes$beta.indexes] = -10^9
    param.upper.bounds[param.indexes$beta.indexes] = 10^9
    if (length(param.indexes$psi.indexes) > 0) {
      # order is age, sex, FACTPal_0
      param.lower.bounds[param.indexes$psi.indexes] = -20
      param.upper.bounds[param.indexes$psi.indexes] = 20
    }
    param.lower.bounds[param.indexes$sigma.index] = 0.00000001
    param.upper.bounds[param.indexes$sigma.index] = 25
    param.lower.bounds[param.indexes$nu.index] = 0.00000001
    param.upper.bounds[param.indexes$nu.index] = 25
    #param.lower.bounds[param.indexes$eta.index] = 0.00000001
    #param.upper.bounds[param.indexes$eta.index] = 25
  }

  # set initial values and bounds for survival model parameters
  if (param.indexes$include.survival) {
    param.lower.bounds[param.indexes$alpha.indexes] = -10^9
    param.upper.bounds[param.indexes$alpha.indexes] = 10^9
    if (length(param.indexes$phi.indexes) > 0) {
      param.lower.bounds[param.indexes$phi.indexes] = -5
      param.upper.bounds[param.indexes$phi.indexes] = 5
    }
  }

  results = list()
  results$param.lower.bounds = param.lower.bounds
  results$param.upper.bounds = param.upper.bounds

  return (results)
}






