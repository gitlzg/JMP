
library(splines)

getSplineBasis = function(knots, # knot values, includes boundary values
    times, # times at which spline should be evaluated
    intercept=TRUE, # include an intercept?
    natural.spline=FALSE # natural or regular cubic spline?
) {
  if (length(knots) == 0) {
    # no spline so return a vector of 1's
    return (as.matrix(rep(1, length(times))))
  } else if (length(knots) == 2) {
    #internal.knots = NULL

    # linear spline return a matrix of 1's and t's
    return (as.matrix(cbind(rep(1, length(times)),times)))

  } else {
    internal.knots = knots[2:(length(knots)-1)]
  }
  boundary.knots = c(knots[1], knots[length(knots)])
  if (natural.spline) {
      basis = ns(x=times, knots=internal.knots, Boundary.knots=boundary.knots, intercept=intercept)
  } else {
      basis = bs(x=times, knots=internal.knots, Boundary.knots=boundary.knots, intercept=intercept, degree=3)
  }
  #if (!intercept) {
  #  basis = cbind(rep(1, nrow(basis)), basis)
  #}

  return (basis)
}

cubicSpline = function(basis, # basis matrix
                       params # basis parameters
) {
    # must cast to matrix in case there is only one basis param
   return (as.matrix(basis) %*% params)
}

computeSplineVar = function(basis, cov.submatrix) {
  #message("Dimensions of basis: ", paste(dim(basis), collapse=","))
  part.1 = basis %*% cov.submatrix
  spline.var = rep(0, nrow(basis))
  for (i in 1:nrow(basis)) {
    spline.var[i] = t(part.1[i,]) %*% basis[i,]
  }
  return (spline.var)
}

getLongitudinalSplineValues = function(long.times, long.knots, param.values, param.indexes, cov.matrix=NA, natural.spline=TRUE) {
  message("Number of non-unique times: ", length(as.vector(as.matrix(long.times))))
  message(paste(sort(as.vector(as.matrix(long.times)))[1:15], collapse=", "))
  unique.times = unique(as.vector(as.matrix(long.times)))
  message("Number of unique times: ", length(unique.times))
  message(paste(sort(unique.times)[1:15], collapse=", "))
  unique.times = unique.times[which(!is.na(unique.times))]
  unique.times = sort(unique.times)
  long.basis = getSplineBasis(knots=long.knots, times=unique.times, natural.spline=natural.spline)
  results = list()
  results$times = unique.times
  results$long.basis = long.basis
  message("ncol long.basis: ", ncol(long.basis), ", length mus: ", length(param.indexes$mu.indexes))
  results$untreated = long.basis %*% param.values[param.indexes$mu.indexes]
  results$treated = results$untreated + long.basis %*% param.values[param.indexes$beta.indexes]
  if (!is.na(cov.matrix)[1]) {
    results$untreated.var = computeSplineVar(long.basis, cov.matrix[param.indexes$mu.indexes, param.indexes$mu.indexes])
    results$treated.var = computeSplineVar(cbind(long.basis, long.basis), cov.matrix[c(param.indexes$mu.indexes, param.indexes$beta.indexes),
            c(param.indexes$mu.indexes, param.indexes$beta.indexes)])
  }
  return (results)
}






