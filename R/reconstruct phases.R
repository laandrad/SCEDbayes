## calculate phase level means from model coefficients

reconstruct.phases <- function(bayes.coeff, model, P) {

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)
  lPhase = c(0,nPhase, nPoints)
  lPhase = lPhase[-1] - lPhase[-length(lPhase)]
  lPhase = lPhase / 2

  N = nrow(bayes.coeff)

  if(model == 'trend'){

      beta1.star = sapply(1:N, function(j) beta1[j] - beta5[j]*lPhase[1] )                                ## B1 level
      beta2.star = sapply(1:N, function(j) beta2[j] - beta6[j]*lPhase[2] )                                ## A2 level
      beta3.star = sapply(1:N, function(j) beta3[j] - beta7[j]*lPhase[3] )                                ## B2 level

      beta0A1 = beta0                                                                                     ## A1 level
      beta0B1 = sapply(1:N, function(j) beta0[j] + beta1.star[j] )                                        ## B1 level
      beta0A2 = sapply(1:N, function(j) beta0[j] + beta1.star[j] + beta2.star[j] )                        ## A2 level
      beta0B2 = sapply(1:N, function(j) beta0[j] + beta1.star[j] + beta2.star[j] + beta3.star[j] )        ## B2 level

      beta1A1 = beta4                                                                                     ## A1 slope
      beta1B1 = sapply(1:N, function(j) beta4[j] + beta5[j] )                                             ## B1 slope
      beta1A2 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j] )                                  ## A2 slope
      beta1B2 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j] + beta7[j] )                       ## B2 slope
      sigma =   sapply(1:N, function(j) sigma.delta[j] )                                                  ## stdev

      return(data.frame(beta0A1, beta0B1, beta0A2, beta0B2, beta1A1, beta1B1, beta1A2, beta1B2, sigma))

  } else{

      beta0A1 = beta0                                                                ## A1 level
      beta0B1 = sapply(1:N, function(j) beta0[j] + beta1[j] )                        ## B1 level
      beta0A2 = sapply(1:N, function(j) beta0[j] + beta1[j]+ beta2[j] )              ## A2 level
      beta0B2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j] + beta3[j] )  ## B2 level
      sigma =   sapply(1:N, function(j) sigma.delta[j] )                             ## stdev

    return(data.frame(beta0A1, beta0B1, beta0A2, beta0B2, sigma))

  }

}



