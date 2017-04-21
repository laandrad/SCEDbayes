## calculate phase level means from model coefficients

reconstruct.coeff <- function(bayes.coeff, model) {

  N = nrow(bayes.coeff)

  if(model == 'trend'){

    beta0 = sapply(1:N, function(j) beta0[j] / (1 - rho[j]) )                     ## A1 level
    beta1 = sapply(1:N, function(j) beta1[j] / (1 - rho[j]) )                     ## B1 level
    beta2 = sapply(1:N, function(j) beta2[j] / (1 - rho[j]) )                     ## A2 level
    beta3 = sapply(1:N, function(j) beta3[j] / (1 - rho[j]) )                     ## B2 level
    beta4 = beta4                                                                 ## A1 slope
    beta5 = beta5                                                                 ## B1 slope
    beta6 = beta6                                                                 ## A2 slope
    beta7 = beta7                                                                 ## B2 slope
    sigma =   sapply(1:N, function(j) sigma.delta[j] / (1 - rho[j]) )             ## stdev

    return(data.frame(beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7, sigma, rho))

  } else{


    beta0 = sapply(1:N, function(j) beta0[j] / (1 - rho[j]) )                     ## A1 level
    beta1 = sapply(1:N, function(j) beta1[j] / (1 - rho[j]) )                     ## B1 level
    beta2 = sapply(1:N, function(j) beta2[j] / (1 - rho[j]) )                     ## A2 level
    beta3 = sapply(1:N, function(j) beta3[j] / (1 - rho[j]) )                     ## B2 level
    sigma =   sapply(1:N, function(j) sigma.delta[j] / (1 - rho[j]) )             ## stdev

    return(data.frame(beta0, beta1, beta2, beta3, sigma, rho))

  }

}
