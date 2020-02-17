

#-------------------------------------------------------------------------------------
# Utility function that returns the indexes of various parameters in the required
# parameter vectors and the names of all parameters.
#
# Includes the following for the longitudinal mixed model:
#
# mu: longitudinal mean value, cubic spline basis functions
# beta: coef of treatment indicator in longitudinal model, cubic spline basis functions
# psi: coef of covariates in longitudinal model, vector
# sigma: random effect variance
# nu: stddev of longitudinal Gaussian process 
# eta: exponential coef in Gaussian process correlation 
#
# Includes the following for the Cox survival model:
#
# alpha: coef of treatment indicator in survival model, cubic spline basis functions
# phi: coef of covariates in survival model, vector
#
# Returns a list with the following elements:
# 
# param.names: Vector with names of all parameters in correct order
#    (mu, beta, psi, sigma, nu, eta, alpha, phi)
# mu.indexes
# beta.indexes
# psi.indexes
# sigma.index
# nu.index
# eta.index
# alpha.indexes
# phi.indexes
#-------------------------------------------------------------------
getParamIndexes = function(
    num.long.fixed.effects, # number of fixed effect covariates
    num.surv.fixed.effects, # number of fixed effect covariates    
    num.long.knots, # number of cubic spline knots in the longitudinal model; if 0, no spline
    num.surv.knots, # number of cubic spline knots in the survival model; if 0, no spline
    natural.spline=F, # T for natural cubic spline, false for standard
    num.events, # number of observed events
    include.survival=T, # true to include the survival parameters
    include.longitudinal=T # true to include the longitudinal parameters
) {  
  
  # Determine the names for cubic spline basis function params in the mixed model for the QoL mean and treatment terms
  # With an intercept, the number of basis functions equals the number of knots (internal plus boundary)
  if (num.long.knots > 0) {
    if (natural.spline) {
      num.long.basis.params = num.long.knots
    } else {
      num.long.basis.params = num.long.knots +2   
    }
    mu.names = c()  
    beta.names = c()
    for (i in 1:num.long.basis.params) {
      mu.names = c(mu.names, paste("mu_", i, sep=""))
      beta.names = c(beta.names, paste("beta_", i, sep=""))
    }
  } else {
    num.long.basis.params = 1
    mu.names = "mu"
    beta.names = "beta"
  }

  # Determine the names for the cubic spline basis function params in the survival model:
  if (num.surv.knots > 0) {    
    if (natural.spline) {
      num.surv.basis.params = num.surv.knots 
    } else {
      num.surv.basis.params = num.surv.knots + 2   
    }
    alpha.names  = c()
    for (i in 1:num.surv.basis.params) {
      alpha.names = c(alpha.names, paste("alpha_", i, sep=""))
    }
  } else {
    num.surv.basis.params = 1
    alpha.names = "alpha"
  }
  
  # Determine the number of covariate coefficients in both the QoL and survival models:
  psi.names = phi.names = c()
  if (num.long.fixed.effects > 0) {
    psi.names = paste("psi", 1:num.long.fixed.effects, sep="")
  } 
  if (num.surv.fixed.effects > 0) {
    phi.names = paste("phi", 1:num.surv.fixed.effects, sep="")  
  }   
  
    
  # Define a vector of all parameter names
  param.names = c()
  if (include.longitudinal) {
    #param.names = c(mu.names, beta.names, psi.names, "sigma", "nu", "eta")
    param.names = c(mu.names, beta.names, psi.names, "sigma", "nu")  
  }  
  if (include.survival) {
    param.names = c(param.names, alpha.names, phi.names)#, lambda.names)
  }
  
  result = list()
  result$include.survival = include.survival
  result$include.longitudinal = include.longitudinal
  result$num.long.basis.params = num.long.basis.params
  result$num.surv.basis.params = num.surv.basis.params
  result$param.names = param.names  
  
  # Define the indices of the different parameter groups in the params vector
  start = 1
  if (include.longitudinal) {
    # mu indexes
    end = num.long.basis.params
    result$mu.indexes = start:end 
    # beta indexes
    start = end+1
    end = end + num.long.basis.params
    result$beta.indexes = start:end
    # psi indexes
    if (num.long.fixed.effects > 0) {
      start = end+1
      end = end + num.long.fixed.effects
      result$psi.indexes = start:end
    } else {
      result$psi.indexes = c() 
    }
    # sigma, tau, nu, eta index
    result$sigma.index = end+1  
    result$nu.index = end+2
    #result$eta.index = end+3
    #start = end+4    
    start = end+3      
  } 
  
  if (include.survival) {    
    # alpha indexes
    end = start + num.surv.basis.params - 1
    result$alpha.indexes = start:end     
    # phi indexes
    if (num.surv.fixed.effects > 0) {
      start = end+1
      end = end + num.surv.fixed.effects
      result$phi.indexes = start:end
    }
  }
  
  return (result)
}



