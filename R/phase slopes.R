## calculate phase level means from model coefficients

phase.trends <- function(bayes.coeff) {

  N = nrow(bayes.coeff)

  beta0A1 = sapply(1:N, function(j) beta0[j])                                   ## A1 level
  beta0B1 = sapply(1:N, function(j) beta0[j] + beta1[j])                        ## B1 level
  beta0A2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j])             ## A2 level
  beta0B2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j] + beta3[j])  ## B2 level
  beta1A1 = beta4                                                               ## A1 slope
  beta1B1 = sapply(1:N, function(j) beta4[j] + beta5[j])                        ## B1 slope
  beta1A2 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j])             ## A2 slope
  beta1B2 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j] + beta7[j])  ## B2 slope
  sigma =   sapply(1:N, function(j) 1/sqrt(tau[j]))                             ## stdev
  
  return(data.frame(beta0A1, beta0B1, beta0A2, beta0B2, beta1A1, beta1B1, beta1A2, beta1B2, sigma))

}
