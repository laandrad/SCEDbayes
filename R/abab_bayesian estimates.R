## Calculate Bayesian estimates
abab.bayesian.estimates <- function(y, P, s, model) {
  phases = abab.create.phases(P, s)

  data = list(
    nPoints = length(y),
    y = y,
    P1 = unlist(phases[1]),
    P2 = unlist(phases[2]),
    P3 = unlist(phases[3]),
    T1 = unlist(phases[4]),
    T2 = unlist(phases[5]),
    T3 = unlist(phases[6]),
    T4 = unlist(phases[7])
  )

  ## Load model
  if(model == 'trend'){
    reg = lm(y ~ T1 + P1 + T2:P1 + P2 + T3:P2 + P3 + T4:P3, data = data)
    initsList = list(
      beta0 = reg$coefficients[1] ,
      beta1 = reg$coefficients[2] ,
      beta2 = reg$coefficients[3] ,
      beta3 = reg$coefficients[4] ,
      beta4 = reg$coefficients[5] ,
      beta5 = reg$coefficients[6] ,
      beta6 = reg$coefficients[7] ,
      beta7 = reg$coefficients[8] ,
      rho = 0,
      nu = 1,
      sigma.delta = length(y) / sum(reg$residuals^2)
    )
    param = c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'beta7', 'rho', 'nu', 'sigma.delta')

  } else {
    reg = lm(y ~ P1 + P2 + P3, data = data)
    initsList = list(
      beta0 = reg$coefficients[1] ,
      beta1 = reg$coefficients[2] ,
      beta2 = reg$coefficients[3] ,
      beta3 = reg$coefficients[4] ,
      rho = 0,
      nu = 1,
      sigma.delta = length(y) / sum(reg$residuals^2)
    )
    param = c('beta0', 'beta1', 'beta2', 'beta3', 'rho', 'nu', 'sigma.delta')
    data = data[1:5]

  }

  ## Initialize MCMC
  nIter = ceiling((numSavedSteps * thinSteps)/nChains)

  jagsModel = jags.model('model.txt', data=data, n.chains = nChains, n.adapt = adaptSteps, inits=initsList)

  ## Burn in chain
  cat('Burning in the MCMC chain...\n')
  update(jagsModel, n.iter=burnInSteps)

  ## Compute final chain
  cat('Sampling final MCMC chain...\n')
  codaSamples = coda.samples(jagsModel, variable.names=param, n.iter=nIter, thin=thinSteps)

  ## Save the sampled chains as a matrix for further processing
  mcmc = as.matrix(codaSamples)

  return(list(mcmc, codaSamples))

}
