## calculate standardized effect sizes for intercept only model

delta.AR1 <- function(y, P, s, bayes.coeff) {

  phases = create.phases(P, s)
  nPoints = length(y)
  P1 = unlist(phases[1])
  P2 = unlist(phases[2])
  P3 = unlist(phases[3])
  N = nrow(bayes.coeff)

  ## calculate the within-subject stdev from residuals = observed - predicted value at time t

  stdev = sapply(1:N, function(j) {
    res = sapply(2:nPoints, function(t){
      yhat = beta0[j]*(1-rho) + beta1[j] * (P1[t]-rho*P1[t-1]) + beta2[j] * (P2[t]-rho*P2[t-1]) + beta3[j] * (P3[t]-rho*P3[t-1])
      yhat = c( beta0[j] + beta1[j] * P1[1] + beta2[j] * P2[1] + beta3[j] * P3[1] , yhat)
      return((y[t] - yhat)^2)
    })
    return(sqrt(sum(res)/nPoints))
  })

  ## calculate standardized effect sizes
  delta1 = sapply(1:N, function(j) beta1[j] / stdev[j]) ## A1B1
  delta2 = sapply(1:N, function(j) beta2[j] / stdev[j]) ## B1A2
  delta3 = sapply(1:N, function(j) beta3[j] / stdev[j]) ## A2B2

  return(data.frame(delta1, delta2, delta3))

}
