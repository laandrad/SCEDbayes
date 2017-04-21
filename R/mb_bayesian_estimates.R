## Calculate Bayesian estimates
mb.bayesian.estimates <- function(y, P, s, model) {
  phases = mb.create.phases(P, s)

  ## Load data
  if(model == 'trend'){
    data = list(
      nPoints = length(y),
      y = y,
      P1 = unlist(phases[1]),
      T1 = unlist(phases[2]),
      T2 = unlist(phases[3])
    )

  } else {

    data = list(
      nPoints = length(y),
      y = y,
      P1 = unlist(phases[1])
    )
  }

  ## Load model
  if(model == 'trend'){
    reg = lm(y ~ T1 + P1 + T2:P1, data = data)
    initsList = list(
      beta0 = reg$coefficients[1] ,
      beta1 = reg$coefficients[2] ,
      beta2 = reg$coefficients[3] ,
      beta3 = reg$coefficients[4] ,
      rho = 0,
      nu = 1,
      sigma.delta = length(y) / sum(reg$residuals^2)
    )
    param = c('beta0', 'beta1', 'beta2', 'beta3', 'rho', 'sigma.delta')

  } else {
    reg = lm(y ~ P1, data = data)
    initsList = list(
      beta0 = reg$coefficients[1] ,
      beta1 = reg$coefficients[2] ,
      rho = 0,
      nu = 1,
      sigma.delta = length(y) / sum(reg$residuals^2)
    )
    param = c('beta0', 'beta1', 'rho', 'nu', 'sigma.delta')

  }

  ## Initialize MCMC
  adaptSteps = 1000
  burnInSteps = 1000
  nChains = 3
  numSavedSteps = 20000
  thinSteps = 1
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
  return(mcmc)

}
