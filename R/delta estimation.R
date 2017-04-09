## calculate standardized effect sizes for intercept only model

compute.delta <- function(y, P, s, bayes.coeff, model) {

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
                beta0[j]*(1-rho[j])             + beta4[j]*(T1[t]       - rho[j]*T1[t-1])         +
                beta1[j]*(P1[t]-rho[j]*P1[t-1]) + beta5[j]*(T2[t]*P1[t] - rho[j]*T2[t-1]*P1[t-1]) +
                beta2[j]*(P2[t]-rho[j]*P2[t-1]) + beta6[j]*(T3[t]*P2[t] - rho[j]*T3[t-1]*P2[t-1]) +
                beta3[j]*(P3[t]-rho[j]*P3[t-1]) + beta7[j]*(T4[t]*P3[t] - rho[j]*T4[t-1]*P3[t-1]) +
                rho[j] * y[t-1])

      yhat = c( beta0[j]         + beta4[j] * T1[1]       +
                beta1[j] * P1[1] + beta5[j] * T2[1]*P1[1] +
                beta2[j] * P2[1] + beta6[j] * T3[1]*P2[1] +
                beta3[j] * P3[1] + beta7[j] * T4[1]*P3[1] ,
                yhat)

      res = (y - yhat)^2

      return(sqrt(sum(res)/nPoints))
    })

    ## calculate effect sizes
    delta1 = sapply(1:N, function(j) beta1[j] - beta5[j]*lPhase[1] )          ## A1B1
    delta2 = sapply(1:N, function(j) beta2[j] - beta6[j]*lPhase[2] )          ## B1A2
    delta3 = sapply(1:N, function(j) beta3[j] - beta7[j]*lPhase[3] )          ## A2B2

    ## standardize effect sizes
    delta1 = sapply(1:N, function(j) delta1[j] / stdev[j])                    ## A1B1
    delta2 = sapply(1:N, function(j) delta2[j] / stdev[j])                    ## B1A2
    delta3 = sapply(1:N, function(j) delta3[j] / stdev[j])                    ## A2B2

    return(data.frame(delta1, delta2, delta3))

  } else{

    ## calculate the within-subject stdev from residuals = observed - predicted value at time t
    stdev = sapply(1:N, function(j) {
      yhat = sapply(2:nPoints, function(t) beta0[j] * (1-rho[j]) +
                                           beta1[j] * (P1[t]-rho[j]*P1[t-1]) +
                                           beta2[j] * (P2[t]-rho[j]*P2[t-1]) +
                                           beta3[j] * (P3[t]-rho[j]*P3[t-1]) )

      yhat = c( beta0[j] + beta1[j] * P1[1] + beta2[j] * P2[1] + beta3[j] * P3[1] , yhat )
      res = (y - yhat)^2
      return(sqrt(sum(res)/nPoints))
    })

    ## calculate standardized effect sizes
    delta1 = sapply(1:N, function(j) beta1[j] / stdev[j]) ## A1B1
    delta2 = sapply(1:N, function(j) beta2[j] / stdev[j]) ## B1A2
    delta3 = sapply(1:N, function(j) beta3[j] / stdev[j]) ## A2B2

    return(data.frame(delta1, delta2, delta3))

  }
}
