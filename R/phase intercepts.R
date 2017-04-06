## calculate phase level means from model coefficients

phase.levels <- function(bayes.coeff) {

  N = nrow(bayes.coeff)

  beta0A1 = beta0                                                               ## A1 level
  beta0B1 = sapply(1:N, function(j) beta0[j] + beta1[j])                        ## B1 level
  beta0A2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j])             ## A2 level
  beta0B2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j] + beta3[j])  ## B2 level
  sigma =   sapply(1:N, function(j) 1/sqrt(tau[j]))                             ## stdev

  return(data.frame(beta0A1, beta0B1, beta0A2, beta0B2, sigma))

}
