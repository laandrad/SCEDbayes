## calculate phase level means from model coefficients

mb.reconstruct.phases <- function(bayes.coeff, model, P) {

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)
  lPhase = c(0,nPhase, nPoints)
  lPhase = lPhase[-1] - lPhase[-length(lPhase)]
  lPhase = (lPhase-1) / 2

  N = nrow(bayes.coeff)

  if(model == 'trend'){

      beta1.star = sapply(1:N, function(j) beta1[j] - beta3[j]*lPhase[1] )                               ## B level

      beta0A = beta0                                                                                     ## A level
      beta0B = sapply(1:N, function(j) beta0[j] + beta1.star[j] )                                        ## B level

      beta1A = beta2                                                                                     ## A slope
      beta1B = sapply(1:N, function(j) beta2[j] + beta3[j] )                                             ## B slope
      sigma =   sapply(1:N, function(j) sigma.delta[j] )                                                 ## stdev

      return(data.frame(beta0A, beta0B, beta1A, beta1B, sigma))

  } else{

      beta0A = beta0                                                                ## A level
      beta0B = sapply(1:N, function(j) beta0[j] + beta1[j] )                        ## B level
      sigma =   sapply(1:N, function(j) sigma.delta[j] )                            ## stdev

    return(data.frame(beta0A, beta0B, sigma))

  }

}



