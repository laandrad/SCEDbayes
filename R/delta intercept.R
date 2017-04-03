## calculate standardized effect sizes for intercept only model

delta.level <- function(y, P, s, bayes.coeff) {

  phases = create.phases(P, s)
  nPoints = length(y)
  P1 = unlist(phases[1])
  P2 = unlist(phases[2])
  P3 = unlist(phases[3])
  N = nrow(bayes.coeff)

  ## calculate the within-subject stdev from residuals = observed - predicted value at time t

  stdev = sapply(1:N, function(j) {
    res = sapply(1:nPoints, function(i){
      yhat = beta0[j] + beta1[j] * P1[i] + beta2[j] * P2[i] + beta3[j] * P3[i]
      return((y[i] - yhat)^2)
    })
    return(sqrt(sum(res)/nPoints))
  })

  ## calculate standardized effect sizes
  delta1 = sapply(1:N, function(j) beta1[j] / stdev[j]) ## A1B1
  delta2 = sapply(1:N, function(j) beta2[j] / stdev[j]) ## B1A2
  delta3 = sapply(1:N, function(j) beta3[j] / stdev[j]) ## A2B2

  return(data.frame(delta1, delta2, delta3))

}
