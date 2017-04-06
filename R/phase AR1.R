## calculate phase level means from model coefficients

phase.AR1 <- function(bayes.coeff) {

  N = nrow(bayes.coeff)

  # sapply(1:N, function(j) print((beta0[j] + beta1[j])/(1-rho[j])))
  beta0A1 = beta0                                                                            ## A1 level
  beta0B1 = sapply(1:N, function(j) (beta0[j] + beta1[j])/(1-rho[j]))                        ## B1 level
  beta0A2 = sapply(1:N, function(j) (beta0[j] + beta1[j]+ beta2[j])/(1-rho[j]))              ## A2 level
  beta0B2 = sapply(1:N, function(j) (beta0[j] + beta1[j] + beta2[j] + beta3[j])/(1-rho[j]))  ## B2 level
  sigma =   sapply(1:N, function(j) (sigma.delta[j]/(1-rho[j])))                                 ## stdev
  
  return(data.frame(beta0A1, beta0B1, beta0A2, beta0B2, sigma))

}



