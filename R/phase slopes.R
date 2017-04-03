## calculate phase level means from model coefficients

phase.trends <- function(bayes.coeff) {

  N = nrow(bayes.coeff)

  gamma0 = sapply(1:N, function(j) beta0[j])                                   ## A1 level
  gamma1 = sapply(1:N, function(j) beta0[j] + beta1[j])                        ## B1 level
  gamma2 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j])             ## A2 level
  gamma3 = sapply(1:N, function(j) beta0[j] + beta1[j] + beta2[j] + beta3[j])  ## B2 level
  gamma4 = beta4                                                               ## A1 slope
  gamma5 = sapply(1:N, function(j) beta4[j] + beta5[j])                        ## B1 slope
  gamma6 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j])             ## A2 slope
  gamma7 = sapply(1:N, function(j) beta4[j] + beta5[j] + beta6[j] + beta7[j])  ## B2 slope

  return(data.frame(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7))

}

