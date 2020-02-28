##' Joint Model for Palliative care studies
##'
##' A semiparametric joint model for terminal trend of quality of life (QOL) and survival
##'
##' - For longitudinal QOL outcomes: semiparametric mixed effects submodel (only random intercept is used in this package)
##'
##'   \deqn{Y_i(t^*) = \beta_\mu(t^*) + A_i \beta_A(t^*) + X_i^T \psi_X + b_i + \epsilon_i(t^*)}
##'
##' - For survival outcomes: Cox submodel
##'   \deqn{\lambda_i(t) = \lambda_0(t) exp(A_i \alpha_A + \tilde{X}_i^T \alpha_X)}
##'
##' @title Joint Model for Palliative care studies
##' @param dat formatted R dataset. See example dataset: data(dat_JMP)
##' @param dataPath path of the data set to read in (if dat not supplied)
##' @param covariate.fields covariates adjusted for
##' @param qol.prefix prefix of the longitudinal quality of life measurements
##' @param qol.time.prefix prefix of the time variable for longitudinal quality of life measurements
##' @param id.field id variable in the dataset
##' @param survival.time.field survival time variable in the dataset
##' @param censoring.status.field censoring variable (1 for death, 0 for censored)
##' @param treatment.status.field treatment status (1 for treated,0 for control)
##' @param long.knot.nums number of knots for nonparametric curves (default 2 to 11)
##' @return A list of analysis results
##'
##'         1. monthly QOL results:
##'
##'           - month: time variable
##'
##'           - nCtr, nTrt: observed number of patients in control and treatment group, respectively
##'
##'           - CtrQOL, TrtQOL: average monthly QOL in control and treatment group
##'
##'           - Ctr95\%CIlow, Ctr95\%CIhi: lower and upper bounds of 95% CI for control group
##'
##'           - Trt95\%CIlow, Trt95\%CIhi: lower and upper bounds of 95% CI for treatment group
##'
##'           - TrtEffect: monthly treatment effect
##'
##'           - pValueTrtEffect: p-value for monthly treatment effect
##'
##'         2. The overall p value for testing the entire model
##'
##'         3. knot.combinations: different combinations of knots for longitudinal and survival submodels with their corresponding AIC and BIC values
##'
##'         4. maximum likelihood estimates (mle) results
##'
##'         and
##'
##'         Four figures saved in the current working directory
##'
##'         1. estimated terminal trend in a png file called "Longitudinal_trajectories.png"
##'
##'         2. the time-varying treatment effect in a png file called "Longitudinal_spline_treatment_effect.png"
##'
##'         3. estimated survival function and in a png file "Survival_function.png"
##'
##'         4. estimated cumulative hazard function in a png file "Cumulative_hazard.png"
##'
##' @author Zhigang Li <zhigang.li@@ufl.edu>
##'
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @references
##' Li Z, Frost HR, Tosteson TD, et al. A Semiparametric Joint Model for Terminal Trend of Quality of Life and Survival in Palliative Care Research. Statistics in Medicine. 2017;36:4692–4704. https://doi.org/10.1002/sim.7445
##' @export
##' @import mefa nlme numDeriv survival splines grDevices graphics stats utils
##' @examples
##' data(dat_JMP)
##'
##' results <- JMP(dat=dat_JMP,
##' covariate.fields = c("sex", "qol_0"),
##' qol.prefix="qol_", qol.time.prefix="time_",
##' id.field = "id", survival.time.field = "survival_time",
##' censoring.status.field = "death", treatment.status.field = "trt")
##'
##' # Extract monthly results from output
##' results$monthlyResults
##'
##' # The overall p value for testing the entire model
##' results$overallPvalueModel


JMP<- function(
  dat=NA,
  dataPath=NA,
  covariate.fields,
  qol.prefix,
  qol.time.prefix,

  id.field,
  survival.time.field,
  censoring.status.field,
  treatment.status.field,

  long.knot.nums=2:11

){



  # 1. JointModeling_ENABLE_Test.R

  #message("Selecting optimal number of knots for ENABLE")

  knot.combinations = suppressMessages(knotSelectionViaBIC(
    optimize=FALSE,
    naive=FALSE,
    include.longitudinal=TRUE,
    include.survival=TRUE,
    long.knot.nums=long.knot.nums,
    surv.knot.nums=0,
    enabledata.path=dataPath,
    enabledata=dat,
    long.covariate.fields = covariate.fields,
    surv.covariate.fields = covariate.fields,
    qol.prefix=qol.prefix,
    qol.time.prefix=qol.time.prefix,
    has.qol.t0 = FALSE,
    reltol=1e-8, # convergence tolerance for optim()
    outer.iterations=10, # number of outer iterations for calls to constrOptim(); done for all methods except "L-BFGS-B"
    outer.eps=1e-8, # convergence tolerance for constrOptim()
    optim.method="BFGS", #Nelder-Mead",
    enforce.bounds=TRUE, # Controls whether contrOptim() or optim() is called
    natural.spline=TRUE,
    sample=NA,
    id.field = id.field,
    survival.time.field=survival.time.field,
    censoring.status.field=censoring.status.field,
    treatment.status.field=treatment.status.field)
)#supressM

  results_JMP <- list(knot.combinations=knot.combinations)

  #print(("AIC/BIC results for sample "))
  ordered.knot.combinations = knot.combinations[order(knot.combinations$AIC),]
  #ordered.knot.combinations = knot.combinations[order(knot.combinations$BIC),]
  num.long.knots = ordered.knot.combinations$num.long.knots[1]
  num.surv.knots = ordered.knot.combinations$num.surv.knots[1]
  #print(paste("Using ", num.long.knots, " longitudinal knots and ", num.surv.knots, " survival knots according to AIC/BIC criterion"))

  results=suppressMessages(testENABLE(
    use.saved.data = FALSE,
    save.data=FALSE,
    data.file = NA,
    optimize=TRUE,
    naive=FALSE,
    compute.standard.errors=TRUE,
    profile=FALSE,
    include.longitudinal=TRUE,
    include.survival=TRUE,
    num.long.knots=num.long.knots,
    num.surv.knots=num.surv.knots,
    enabledata.path=dataPath,
    long.covariate.fields = covariate.fields,
    surv.covariate.fields = covariate.fields,
    qol.prefix=qol.prefix,
    qol.time.prefix=qol.time.prefix,
    has.qol.t0 = FALSE,
    reltol=1e-8, # convergence tolerance for optim()
    outer.iterations=10, # number of outer iterations for calls to constrOptim(); done for all methods except "L-BFGS-B"
    outer.eps=1e-8, # convergence tolerance for constrOptim()
    optim.method="BFGS", #Nelder-Mead",
    enforce.bounds=TRUE, # Controls whether contrOptim() or optim() is called
    natural.spline=TRUE,
    sample=NA,
    resample.method=NA,
    enabledata.table=dat,
    id.field = id.field,
    survival.time.field=survival.time.field,
    censoring.status.field=censoring.status.field,
    treatment.status.field=treatment.status.field)
)#supressM

  # save these so no need to output rData file
   mle_save=results$mle
   init.value.results_save=results$init.value.results
   enable.data_save=results$enable.data


  # 2. JointModeling_ENABLE_Analyze.R


  #----------------------------------------------------------------------------
  # Analysis parameters
  #----------------------------------------------------------------------------


  long.spline.file = "Longitudinal_trajectories.png"
  treatment.effect.spline.file = "Longitudinal_spline_treatment_effect.png"
  surv.spline.file = "Survival_spline.png"
  surv.file = "Survival_function.png"
  cum.hazard.file = "Cumulative_hazard.png"
  figure.units="in"
  figure.width=7
  figure.height=5
  figure.res=500

  monthTimePoints=8
  multiTestTimePoints=c(1,2,3)

  #----------------------------------------------------------------------------
  # Define a utility function for computing the longitudinal splines
  #----------------------------------------------------------------------------

  getLongitudinalSplineValues = function(long.times,
                                         avg.long.covariates, #.treated, avg.long.covariates.untreated,
                                         long.knots, param.values, param.indexes, cov.matrix=NA, natural.spline=TRUE) {
    message("Number of non-unique times: ", length(as.vector(as.matrix(long.times))))
    message(paste(sort(as.vector(as.matrix(long.times))), collapse=", "))
    unique.times = unique(as.vector(as.matrix(long.times)))
    message("Number of unique times: ", length(unique.times))
    message(paste(sort(unique.times), collapse=", "))
    unique.times = unique.times[which(!is.na(unique.times))]
    unique.times = sort(unique.times)
    long.basis = getSplineBasis(knots=long.knots, times=unique.times, natural.spline=natural.spline)
    results = list()
    results$times = unique.times
    results$long.basis = long.basis
    message("num long knots: ", length(long.knots))
    message("ncol long.basis: ", ncol(long.basis), ", length mus: ", length(param.indexes$mu.indexes))
    # Compute just the nonparametric component at each time point for treated and untreated
    results$untreated = long.basis %*% param.values[param.indexes$mu.indexes]
    results$treat.effect = long.basis %*% param.values[param.indexes$beta.indexes]
    results$treated = results$untreated + results$treat.effect
    # Compute the average of the parameteric components for treated and untreated
    avg.cov = t(avg.long.covariates) %*% param.values[param.indexes$psi.indexes]
    message("Average covariate value: ", avg.cov)
    #avg.cov.untreated = t(avg.long.covariates.untreated) %*% param.values[param.indexes$psi.indexes]
    #message("Average covariate value for untreated: ", avg.cov.untreated)
    #avg.cov.treated = t(avg.long.covariates.treated) %*% param.values[param.indexes$psi.indexes]
    #message("Average covariate value for treated: ", avg.cov.treated)
    # Add the parametric to the nonparametric
    results$untreated = results$untreated + rep(avg.cov, length(results$untreated))
    results$treated = results$treated + rep(avg.cov, length(results$treated))
    if (!is.na(cov.matrix)[1]) {
      covariate.mat = matrix(rep(avg.long.covariates, nrow(long.basis)), nrow=nrow(long.basis), byrow=TRUE)
      untreated.indexes = c(param.indexes$mu.indexes, param.indexes$psi.indexes)
      results$untreated.var = computeSplineVar(cbind(long.basis, covariate.mat), cov.matrix[untreated.indexes, untreated.indexes])
      treated.indexes = c(param.indexes$mu.indexes, param.indexes$beta.indexes, param.indexes$psi.indexes)
      results$treated.var = computeSplineVar(cbind(long.basis, long.basis, covariate.mat), cov.matrix[treated.indexes, treated.indexes])
      treat.effect.indexes = c(param.indexes$beta.indexes)
      results$treat.effect.var = computeSplineVar(long.basis, cov.matrix[treat.effect.indexes, treat.effect.indexes])
      diff.indexes = c(param.indexes$beta.indexes, param.indexes$psi.indexes)
      results$diff.var = computeSplineVar(cbind(long.basis, covariate.mat), cov.matrix[diff.indexes, diff.indexes])
    }
    return (results)
  }

  #
  ## terminal trend for delayed start group in the 2013 Stats in Med paper
  #
  # QoL.untreated.StatsInMed=function(t){
  #   mean.t=0.6317*126.1855+23.0708+6*(4.214-0.0044)*(t>6)+4.214*t*(t<=6)+0.0044*t*(t>6)
  #   return(mean.t)
  # }

  # terminal trend for early start group
  # QoL.treated.StatsInMed=function(t){
  #   mean.c=0.6317*126.1855+23.0708+6*(4.214-0.0044)*(t>6)+4.214*t*(t<=6)+0.0044*t*(t>6)+
  #       5.8569+6*(-0.6649-(-0.0357))*(t>6)+(-0.6649)*t*(t<=6)+(-0.0357)*t*(t>6)
  #   return(mean.c)
  # }


  #----------------------------------------------------------------------------
  # Load the previously saved estimation results
  #----------------------------------------------------------------------------



  results=suppressMessages(testENABLE(
    use.saved.data = TRUE,
    save.data=TRUE,
    data.file=NA,
    optimize=TRUE,
    naive=FALSE,
    compute.standard.errors=TRUE,
    profile=FALSE,
    include.longitudinal=TRUE,
    include.survival=TRUE,
    num.long.knots=NA, # this is automatically loaded from saved results
    num.surv.knots=NA,# this is automatically loaded from saved results
    enabledata.path=dataPath,
    long.covariate.fields = covariate.fields,
    surv.covariate.fields = covariate.fields,
    qol.prefix=qol.prefix,
    qol.time.prefix=qol.time.prefix,
    has.qol.t0 = FALSE,
    reltol=1e-8, # convergence tolerance for optim()
    outer.iterations=10, # number of outer iterations for calls to constrOptim(); done for all methods except "L-BFGS-B"
    outer.eps=1e-8, # convergence tolerance for constrOptim()
    optim.method="BFGS", #Nelder-Mead",
    enforce.bounds=TRUE, # Controls whether contrOptim() or optim() is called
    natural.spline=TRUE,
    sample=NA,
    resample.method=NA,enabledata.table=dat,
    mle=mle_save,
    init.value.results=init.value.results_save,
    enable.data_load=enable.data_save,

    id.field = id.field,
    survival.time.field=survival.time.field,
    censoring.status.field=censoring.status.field,
    treatment.status.field=treatment.status.field)
)#supressM
  #----------------------------------------------------------------------------
  # Get the parameter covariance matrix
  #----------------------------------------------------------------------------
  cov.matrix = results$mle$cov.mat

  #----------------------------------------------------------------------------
  # Print the results
  #----------------------------------------------------------------------------

  #print(results$mle)
  #print(results$long.knots)
  #print(results$surv.knots)

  results_JMP$mle <- results$mle

  #
  ## compute the covariate average for treated and untreated
  #

  avg.long.covariates = apply(results$enable.data$long.covariates, 2, mean)
  #avg.long.covariates.treated = apply(results$enable.data$long.covariates[which(results$enable.data$treatment.status == 1),], 2, mean)
  #avg.long.covariates.untreated = apply(results$enable.data$long.covariates[which(results$enable.data$treatment.status == 0),], 2, mean)


  #----------------------------------------------------------------------------
  # Plot the longitudinal spline (untreated and treated)
  #----------------------------------------------------------------------------

  png(filename=long.spline.file, width=figure.width, height=figure.height, units=figure.units, res=figure.res)

  spline.values =suppressMessages( getLongitudinalSplineValues(long.times=results$enable.data$long.times,
                                              avg.long.covariates = avg.long.covariates,
                                              long.knots=results$long.knots,
                                              param.values=results$mle$optim.results$par,
                                              param.indexes=results$param.indexes,
                                              cov.matrix=cov.matrix)
  )#supressM

  # piecewise.linear.untreated = QoL.untreated.StatsInMed(spline.values$times)
  # piecewise.linear.treated = QoL.treated.StatsInMed(spline.values$times)
  upper.ci.untreated = spline.values$untreated + 1.96*sqrt(spline.values$untreated.var)
  lower.ci.untreated = spline.values$untreated - 1.96*sqrt(spline.values$untreated.var)
  upper.ci.treated = spline.values$treated + 1.96*sqrt(spline.values$treated.var)
  lower.ci.treated = spline.values$treated - 1.96*sqrt(spline.values$treated.var)
  upper.ci.diff = spline.values$treated + 1.96*sqrt(spline.values$treat.effect.var)
  lower.ci.diff = spline.values$treated - 1.96*sqrt(spline.values$treat.effect.var)
  all.values = c(spline.values$treated, spline.values$untreated,
                 upper.ci.untreated,lower.ci.untreated,
                 upper.ci.treated,lower.ci.treated,
                 upper.ci.diff, lower.ci.diff)
  ylims = c(min(all.values), max(all.values))
  xlims = c(min(c(spline.values$times,results$long.knots)), max(c(spline.values$times, results$long.knots)))
  par(mar=c(4,4,1,1))
  plot(spline.values$times, spline.values$untreated, ylim=ylims, xlim=xlims,
       cex=1, cex.main=.8, cex.lab=.7, cex.axis=.7,
       #main="Estimated QOL scores",
       xlab="QOL measurement time (months retrospective from death)",

       ylab="QOL measurements",
       type="l", col="black", lty="dashed", pch=20)
  #message("Length of untreated: ", length(spline.values$untreated), ", length var: ", length(spline.values$untreated.var))
  lines(spline.values$times, spline.values$treated, col="black", pch=20, lty="solid")
  lines(spline.values$times, upper.ci.diff, col="black", lty="dotted")
  lines(spline.values$times, lower.ci.diff, col="black", lty="dotted")


  legend("bottomright", inset=c(0.03,0.03),
         c("treatment", "control",
           "95% CI of treatment difference centered on treatment group"),
         col=c("black", "black", "black"), lty=c("solid", "dashed", "dotted"), cex=.6)

  dev.off()

  #----------------------------------------------------------------------------
  # Plot the longitudinal spline (treatment effect)
  #----------------------------------------------------------------------------

  png(filename=treatment.effect.spline.file, width=figure.width, height=figure.height, units=figure.units, res=figure.res)

  upper.ci.treat.effect = spline.values$treat.effect + 1.96*sqrt(spline.values$treat.effect.var)
  lower.ci.treat.effect = spline.values$treat.effect - 1.96*sqrt(spline.values$treat.effect.var)
  avg.treat.effect = mean(spline.values$treat.effect)
  all.values = c(spline.values$treat.effect,
                 upper.ci.treat.effect, lower.ci.treat.effect)
  ylims = c(min(all.values), max(all.values))
  xlims = c(min(c(spline.values$times,results$long.knots)), max(c(spline.values$times, results$long.knots)))
  plot(spline.values$times, spline.values$treat.effect, ylim=ylims, xlim=xlims,
       cex=1, cex.main=.8, cex.lab=.7, cex.axis=.7,
       #main=paste("Estimated treatment effect on QOL scores (avg. ", format(avg.treat.effect, digits=3), ")",sep=""),
       xlab="QOL measurement time (months retrospective from death)",
       ylab="QOL measurements",
       type="l", col="black", lty="solid", pch=20)
  abline(h=0, col="black", lty="dashed")
  lines(spline.values$times, upper.ci.treat.effect, col="black", lty="dotted")
  lines(spline.values$times, lower.ci.treat.effect, col="black", lty="dotted")
  legend("bottomright", c(
    expression(paste("Treatment effect (",hat(beta)[A](t*"*"), ",", sep="")),
    expression("95% CI of treatment effect")),
    col=c("black", "black"), lty=c("solid", "dotted"), cex=.6)

  dev.off()

  #----------------------------------------------------------------------------
  # Plot the alpha spline
  #----------------------------------------------------------------------------

  #png(filename=surv.spline.file, width=figure.width, height=figure.height, units=figure.units, res=figure.res)
  ##times = sort(unique(results$enable.data$surv.times))
  #times = seq(from=0, to=18, by=.1)
  #surv.basis = getSplineBasis(knots=results$surv.knots, times=times, natural.spline=T)
  #spline.estimates = results$mle$optim.results$par[results$param.indexes$alpha.indexes]
  #spline.var = computeSplineVar(surv.basis, cov.matrix[results$param.indexes$alpha.indexes, results$param.indexes$alpha.indexes])
  ##message("Alpha indexes: ", paste(results$param.indexes$alpha.indexes, collapse=","))
  ##message(paste(spline.estimates, collapse=","))
  #alpha.values = surv.basis %*% spline.estimates
  #upper.ci.alpha = alpha.values + 1.96*sqrt(spline.var)
  #lower.ci.alpha = alpha.values - 1.96*sqrt(spline.var)
  #all.values = c(alpha.values, upper.ci.alpha, lower.ci.alpha)
  #ylims = c(min(all.values), max(all.values))
  #xlims = c(min(times), max(times))
  #plot(times, alpha.values, ylim=ylims, xlim=xlims, cex.main=.9, cex.lab=.9,
  #    main="Difference between early and delayed treatment on log hazard ratio",
  #    xlab="Time from enrollment (months)",
  #    ylab="Log hazard ratio", type="l", col="black", lty="solid", pch=20)
  #lines(times, upper.ci.alpha, lty="dotted")
  #lines(times, lower.ci.alpha, lty="dotted")
  ##internal.knots = results$long.knots[2:(length(results$long.knots)-1)]
  ##boundary.knots = c(results$long.knots[1], results$long.knots[length(results$long.knots)])
  #abline(h=0, lty="dashed")
  ##abline(v=boundary.knots, lty="dotted")
  ##legend("bottomright", c("Early", "Late"), col=c("black", "black"), lty=c("dotted", "solid"), cex=.8, title="Intervention time")
  #legend("bottomright", c("Early", "Late"), col=c("black", "blue"), lty=c("solid", "solid"), cex=.8, title="Intervention time")
  #
  #dev.off()

  #----------------------------------------------------------------------------
  # Compute cumulative hazard and survival function for treated and untreated
  #----------------------------------------------------------------------------

  # Compute the survival basis for all original survival times
  #surv.basis = getSplineBasis(knots=results$surv.knots, times=results$enable.data$surv.times, natural.spline=T)
  surv.basis=NA
  #message("Dimensions of surv.basis: ", paste(dim(surv.basis), collapse=","))

  # Compute the baseline hazard values
  #message("Computing baseline hazard steps...")
  baseline = suppressMessages(updateBaselineHazardSteps(censored=results$enable.data$censored,
                                       surv.times=results$enable.data$surv.times,
                                       treatment.status=results$enable.data$treatment.status,
                                       surv.basis=surv.basis,
                                       covariates=results$enable.data$surv.covariates,
                                       alphas=results$mle$optim.results$par[results$param.indexes$alpha.indexes],
                                       phis=results$mle$optim.results$par[results$param.indexes$phi.indexes])
  )#supressM
  #message("Finished computing baseline hazard steps.")

  #message("Length of baseline: ", length(baseline))
  #message("Baseline hazard: ", paste(baseline, collapse=","))

  # Get unique uncensored survival times and sort
  uncensored.surv.times = results$enable.data$surv.times[which(results$enable.data$censored ==0)]
  unique.surv.times = unique(uncensored.surv.times)
  #message("Number of unique uncensored surv times: ", length(unique.surv.times))
  sorted.unique.surv.times = sort(unique.surv.times)

  # Sort and subset the survival basis for just the unique survival times
  #basis.indices.to.keep = unique(sapply(sorted.unique.surv.times, function(x) {which(results$enable.data$surv.times == x)[1]}))
  #message("Basis indices to keep: ", paste(basis.indices.to.keep, collapse=","))
  #surv.basis = surv.basis[basis.indices.to.keep,]
  #surv.basis = as.matrix(surv.basis)
  #message("Dimensions of surv.basis: ", paste(dim(surv.basis), collapse=","))
  #message("Dimensions of sorted and filtered surv.basis: ", paste(dim(surv.basis), collapse=","))

  # Compute mean covariate values
  mean.long.cov = apply(results$enable.data$long.covariates, 2, mean)
  mean.surv.cov = apply(results$enable.data$surv.covariates, 2, mean)

  #message("Alpha indexes: ", paste(results$param.indexes$alpha.indexes, collapse=","))

  # Compute the cumulative hazard for treated status
  treated.cum.hazard = rep(0, length(sorted.unique.surv.times))
  untreated.cum.hazard = rep(0, length(sorted.unique.surv.times))
  for (i in 1:length(sorted.unique.surv.times)) {
    surv.time = sorted.unique.surv.times[i]
    treated.hazard.exp = suppressMessages(computeHazardExp(treatment.status=1,
                                          alphas=results$mle$optim.results$par[results$param.indexes$alpha.indexes],
                                          phis=results$mle$optim.results$par[results$param.indexes$phi.indexes],
                                          covariates=mean.surv.cov,
                                          #surv.basis=surv.basis[i,])
                                          surv.basis=NA)
    )#suppress
    #message("Alphas: ", paste(results$mle$optim.results$par[results$param.indexes$alpha.indexes], collapse=","))
    #message("Phis: ", paste(results$mle$optim.results$par[results$param.indexes$phi.indexes], collapse=","))
    #message("Covs: ", paste(mean.surv.cov, collaspe=","))
    untreated.hazard.exp = suppressMessages(computeHazardExp(treatment.status=0,
                                            alphas=results$mle$optim.results$par[results$param.indexes$alpha.indexes],
                                            phis=results$mle$optim.results$par[results$param.indexes$phi.indexes],
                                            covariates=mean.surv.cov,
                                            #surv.basis=surv.basis[i,])
                                            surv.basis=NA)
    )#suppress
    #message("Treated hazard exp for surv time: ", surv.time, ": ", treated.hazard.exp)
    #message("Untreated hazard exp for surv time: ", surv.time, ": ", untreated.hazard.exp)
    if (i == 1) {
      prior.treated.cum.hazard = NA
      prior.untreated.cum.hazard = NA
    } else {
      prior.treated.cum.hazard = treated.cum.hazard[i-1]
      prior.untreated.cum.hazard = untreated.cum.hazard[i-1]
    }
    treated.cum.hazard[i] = suppressMessages(computeCumulativeHazard(surv.time = surv.time,
                                                    hazard.exp.val=exp(treated.hazard.exp),
                                                    event.times=sorted.unique.surv.times,
                                                    baseline.hazard.steps=baseline,
                                                    prior.cum.hazard=prior.treated.cum.hazard)
    )#suppressM
    untreated.cum.hazard[i] = suppressMessages(computeCumulativeHazard(surv.time = surv.time,
                                                      hazard.exp.val=exp(untreated.hazard.exp),
                                                      event.times=sorted.unique.surv.times,
                                                      baseline.hazard.steps=baseline,
                                                      prior.cum.hazard=prior.untreated.cum.hazard)
    )#suppress
    #message("Treated cumulative hazard for surv time: ", surv.time, ": ", treated.cum.hazard[i])
  }

  # Compute the survival function
  treated.surv.func = exp(-treated.cum.hazard)
  untreated.surv.func = exp(-untreated.cum.hazard)
  #message("Length of treated.surv.func: ", length(treated.surv.func))

  #----------------------------------------------------------------------------
  # Plot survival function
  #----------------------------------------------------------------------------

  png(filename=surv.file, width=figure.width, height=figure.height, units=figure.units, res=figure.res)
  plot(sorted.unique.surv.times, treated.surv.func, type="s", pch=20, main="Estimated survival functions", xlab="Time from enrollment (months)", ylab="Survival probability",
       cex=1, cex.main=.8, cex.lab=.7, cex.axis=.7)
  lines(sorted.unique.surv.times, untreated.surv.func, type="s", pch=20, lty="dotted")
  legend("topright", c("treatment", "control"), col=c("black", "black"), lty=c("solid", "dotted"), cex=.6)# title="Intervention time")
  dev.off()

  #----------------------------------------------------------------------------
  # Plot cumulative hazard
  #----------------------------------------------------------------------------

  png(filename=cum.hazard.file, width=figure.width, height=figure.height, units=figure.units, res=figure.res)
  plot(sorted.unique.surv.times, treated.cum.hazard, type="s", pch=20, main="Estimated cumulative hazard", xlab="Time from enrollment (months)", ylab="Cumulative hazard",
       cex=1, cex.main=.8, cex.lab=.7, cex.axis=.7)
  lines(sorted.unique.surv.times, untreated.cum.hazard, type="s", pch=20, lty="dotted")
  legend("bottomright", c("treatment", "control"), col=c("black", "black"), lty=c("solid", "dotted"), cex=.6)# title="Intervention time")
  dev.off()


  #----------------------------------------------------------------------------
  ## to get the values at a vector of retrospective monthly time points
  #----------------------------------------------------------------------------

  # calculate splines values and CI's
  monthlyTime=seq(1,monthTimePoints)
  monthlyValue=suppressMessages(getLongitudinalSplineValues(long.times=monthlyTime,
                                           avg.long.covariates= avg.long.covariates, #.treated, avg.long.covariates.untreated,
                                           long.knots=results$long.knots,
                                           param.values=results$mle$optim.results$par,
                                           param.indexes=results$param.indexes,
                                           cov.matrix=cov.matrix,
                                           natural.spline=TRUE)
  )#suppress

  low.CI.untreat=monthlyValue$untreated-1.96*sqrt(monthlyValue$untreated.var)
  hi.CI.untreat=monthlyValue$untreated+1.96*sqrt(monthlyValue$untreated.var)
  low.CI.treat=monthlyValue$treated-1.96*sqrt(monthlyValue$treated.var)
  hi.CI.treat=monthlyValue$treated+1.96*sqrt(monthlyValue$treated.var)
  pValue.treatEffect=2*(1-pnorm(abs(monthlyValue$treat.effect/sqrt(monthlyValue$treat.effect.var))))

  # calculate sample sizes at each time point
  monthsNo=length(monthlyTime)
  n=length(results$enable.data$surv.times)

  tMatrix=matrix(rep(monthlyTime,n),nrow=n,byrow=TRUE)
  survTMatrix=matrix(rep(results$enable.data$surv.times,monthsNo),nrow=n,byrow=FALSE)
  indicatorM=1*(survTMatrix>=tMatrix)

  nByMonth=matrix(rep(c(NA,NA),monthsNo),nrow=monthsNo)
  for (i in 1:monthsNo) {
    tablei=table(results$enable.data$treatment.status,indicatorM[,i])
    #print(tablei)
    nByMonth[i,1]=tablei[1,2]
    nByMonth[i,2]=tablei[2,2]
  }

  # combine values and sample sizes
  monthlyResults=cbind(monthlyTime,nByMonth,monthlyValue$untreated,low.CI.untreat,hi.CI.untreat,monthlyValue$treated,low.CI.treat,hi.CI.treat,
                       monthlyValue$treat.effect,pValue.treatEffect)

  colnames(monthlyResults)=c("month","nCtr","nTrt","CtrQOL","Ctr95%CIlow","Ctr95%CIhi",
                             "TrtQOL","Trt95%CIlow","Trt95%CIhi",
                             "TrtEffect","pValueTrtEffect")

  #----------------------------------------------------------------------------
  ## to get an overall p value for testing the model and multiple monthly time points
  #----------------------------------------------------------------------------

  #
  ## overall p value for testing the model
  #
  # calculate the vector effect size
  param.values=results$mle$optim.results$par
  betaVector=param.values[results$param.indexes$beta.indexes]

  # calculate the variance of the vector effect
  treat.effect.indexes = c(results$param.indexes$beta.indexes)
  betaSigma=cov.matrix[treat.effect.indexes, treat.effect.indexes]

  chiSquareBeta=t(betaVector)%*%solve(betaSigma)%*%betaVector
  overallPvalueModel=round(1-pchisq(chiSquareBeta, df=length(results$long.knots)),4)


  #
  ## overall p value for testing  multiple monthly time points
  #

  # condition check

  if(length(multiTestTimePoints)<length(results$long.knots)){
    # to calculate the basis at the multiple time points
    multiTimePoints=multiTestTimePoints

    multiTimePoinValues=suppressMessages(getLongitudinalSplineValues(long.times=multiTimePoints,
                                                    avg.long.covariates= avg.long.covariates, #.treated, avg.long.covariates.untreated,
                                                    long.knots=results$long.knots,
                                                    param.values=results$mle$optim.results$par,
                                                    param.indexes=results$param.indexes,
                                                    cov.matrix=cov.matrix,
                                                    natural.spline=TRUE)
    )#suppress

    # calculate chisquare test statistic
    #note: the number of time points cannot be bigger than the number of basis, otherwise,multiSigma is not invertable

    multiBeta=(multiTimePoinValues$long.basis)%*%betaVector
    multiSigma=(multiTimePoinValues$long.basis)%*%betaSigma%*%t(multiTimePoinValues$long.basis)
    chiSquareMulti=t(multiBeta)%*%solve(multiSigma)%*%multiBeta

    overallPvalueMulti=round(1-pchisq(chiSquareMulti, df=length(multiTimePoints)),4)

    # message("The overall p value is ",overallPvalueMulti," for simutaneously testing ", paste(multiTimePoints,","),
    #        " months with a Chi-square test with degrees of freedom ",length(multiTestTimePoints))

  }

  #
  #message("The overall p value is ",overallPvalueModel," for testing the entire model")

  #print(monthlyResults)

  results_JMP$overallPvalueModel <- overallPvalueModel  # 3rd
  results_JMP$monthlyResults <- monthlyResults          # 3th element
  results_JMP <- results_JMP[c(4,3,1,2)]
  return(results_JMP)
}
