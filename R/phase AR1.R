## calculate phase level means from model coefficients

phase.AR1 <- function(bayes.coeff) {

  N = nrow(bayes.coeff)

  # sapply(1:N, function(j) print((beta0[j] + beta1[j])/(1-rho[j])))
  gamma0 = beta0                                                                            ## A1 level
  gamma1 = sapply(1:N, function(j) (beta0[j] + beta1[j])/(1-rho[j]))                        ## B1 level
  gamma2 = sapply(1:N, function(j) (beta0[j] + beta1[j]+ beta2[j])/(1-rho[j]))              ## A2 level
  gamma3 = sapply(1:N, function(j) (beta0[j] + beta1[j] + beta2[j] + beta3[j])/(1-rho[j]))  ## B2 level

  return(data.frame(gamma0, gamma1, gamma2, gamma3))

}



