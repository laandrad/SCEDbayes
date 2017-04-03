## calculate standardized effect sizes for intercept and slopes model

delta.trend <- function(y, P, s, bayes.coeff) {

  phases = create.phases(P, s)
  nPoints = length(y)
  P1 = unlist(phases[1])
  P2 = unlist(phases[2])
  P3 = unlist(phases[3])
  T1 = unlist(phases[4])
  T2 = unlist(phases[5])
  T3 = unlist(phases[6])
  T4 = unlist(phases[7])
  N = nrow(bayes.coeff)

  ## calculate the within-subject stdev from residuals = observed - predicted value at time t
  stdev = sapply(1:N, function(j) {
    res = sapply(1:nPoints, function(i){
      yhat = beta0[j] + beta4[j] * T1[i] + beta1[j] * P1[i] + beta5[j] * T2[i] * P1[i] +
        beta2[j] * P2[i] + beta6[j] * T3[i] * P2[i] + beta3[j] * P3[i] + beta7[j] * T4[i] * P3[i]
      return((y[i] - yhat)^2)
    })
    return(sqrt(sum(res)/nPoints))
  })

  ## calculate effect sizes
  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)
  nPhase = c(nPhase, nPoints)
  nPhase = nPhase[-1] - nPhase[-4]

  delta1 = sapply(1:N, function(j) beta1[j] + (beta5[j] * (nPhase[1]-1))/2 ) ## A1B1
  delta2 = sapply(1:N, function(j) beta2[j] + (beta6[j] * (nPhase[2]-1))/2 ) ## B1A2
  delta3 = sapply(1:N, function(j) beta3[j] + (beta7[j] * (nPhase[3]-1))/2 ) ## A2B2

  ## standardize effect sizes
  delta1 = sapply(1:N, function(j) delta1[j] / stdev[j]) ## A1B1
  delta2 = sapply(1:N, function(j) delta2[j] / stdev[j]) ## B1A2
  delta3 = sapply(1:N, function(j) delta3[j] / stdev[j]) ## A2B2

  return(data.frame(delta1, delta2, delta3))

}
