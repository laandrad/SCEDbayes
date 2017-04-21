## calculate standardized effect sizes for intercept only model

mb.compute.delta <- function(y, P, s, bayes.coeff, model) {

  phases = mb.create.phases(P, s)
  nPoints = length(y)
  P1 = unlist(phases[1])
  T1 = unlist(phases[2])
  T2 = unlist(phases[3])
  N = nrow(bayes.coeff)

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)
  lPhase = c(0,nPhase, nPoints)
  lPhase = lPhase[-1] - lPhase[-length(lPhase)]
  lPhase = (lPhase-1) / 2

  if(model == 'trend'){

    ## calculate the within-subject stdev from residuals = observed - predicted value at time t
    stdev = sapply(1:N, function(j) {

      yhat = sapply(2:nPoints, function(t)
                beta0[j]*(1-rho[j])             + beta2[j]*(T1[t]       - rho[j]*T1[t-1])         +
                beta1[j]*(P1[t]-rho[j]*P1[t-1]) + beta3[j]*(T2[t]*P1[t] - rho[j]*T2[t-1]*P1[t-1]) +
                rho[j] * y[t-1])

      yhat = c( beta0[j]         + beta2[j] * T1[1]       +
                beta1[j] * P1[1] + beta3[j] * T2[1]*P1[1] ,
                yhat)

      res = (y - yhat)^2

      return(sqrt(sum(res)/nPoints))
    })

    ## calculate effect sizes
    delta = sapply(1:N, function(j) beta1[j] - beta3[j]*lPhase[1] )          ## AB

    ## standardize effect sizes
    delta = sapply(1:N, function(j) delta[j] / stdev[j])                    ## AB

    return(delta)

  } else{

    ## calculate the within-subject stdev from residuals = observed - predicted value at time t
    stdev = sapply(1:N, function(j) {
      yhat = sapply(2:nPoints, function(t) beta0[j]*(1-rho[j]) +
                                           beta1[j] * (P1[t]-rho[j]*P1[t-1]) )

      yhat = c( beta0[j] + beta1[j] * P1[1] , yhat )
      res = (y - yhat)^2
      return(sqrt(sum(res)/nPoints))
    })

    ## calculate standardized effect sizes
    delta = sapply(1:N, function(j) beta1[j] / stdev[j]) ## AB

    return(delta)

  }
}
