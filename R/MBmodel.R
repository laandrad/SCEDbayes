## Multiple-Baseline Design

MBmodel = function(y, P, s, model = 'level', plots = TRUE, diagnostics = FALSE, adaptSteps = 5000,
                   burnInSteps = 50000, nChains = 3, numSavedSteps = 100000, thinSteps = 5) {

  ## load packages
  if(!require(rjags)){
    install.packages("rjags")
    library(rjags)
  }

  if(!require(coda)){
    install.packages("coda")
    library(coda)
  }

  ## load model as specified in model argument
  if(model == 'trend'){
    ITS = mb.model.trend
    } else{
    ITS = mb.model.level
  }
  writeLines(ITS, con='model.txt')

  ## calculate bayesian coefficients using a Montecarlo Markov Chain
  beta = mb.bayesian.estimates(y, P, s, model, adaptSteps, burnInSteps, nChains, numSavedSteps, thinSteps)
  chains = beta[[2]]
  beta = data.frame(beta[[1]])

  # return(beta)

  attach(beta)

  ## calculate phase intercept and slope at each point in the mcmc
  cat('Computing phase parameters...\n')

  gamma = mb.reconstruct.phases(beta, model, P)

  cat("  |**************************************************| 100%\n")

  ## calculate effect size for each phase change according to model
  cat('Computing effect size...\n')

  delta = mb.compute.delta(y, P, s, beta, model)

  cat("  |**************************************************| 100%\n")

  detach(beta)

  ## Plotting results
  graphics.off()

  if(plots == TRUE){
    cat('Plotting results...\n')

    plot.titles = c('Standardized Effect Size')

    parameter = c('delta')

    openGraph(width = 14, height = 4)
    layout(matrix(c(1,rep(2,2)), nrow = 1, byrow=T))
    sapply(1, function(i) posterior.plot(delta, plot.titles[i], parameter[i]))

    if(model == 'trend'){
      mb.its.plot.trend(y, P, s, gamma)
    } else{
      mb.its.plot(y, P, s, gamma)
    }

    cat("  |**************************************************| 100%\n")

  } else{
    cat('Plots omitted...\n')
  }

  ## Printing results
  beta.results = t(sapply(beta, describe))
  gamma.results = t(sapply(gamma, describe))
  delta.results = t(data.frame(describe(delta)))
  row.names(delta.results)[1] = 'delta'

  cat('\nBayesian estimates for A and B-A phase change:\n')
  print(beta.results)
  cat('\nRegression estimates for A and B phases:\n')
  print(gamma.results)
  cat('\nStandardized effect size estimates for B-A phase change:\n')
  print(delta.results)

  if(diagnostics == T){
    ## calculate diagnostic statistics
    cat('\nComputing diagnostic statistics...\n')

    openGraph(width = 8, height = 7)
    layout(1)
    gelman.plot(chains)
    GB = gelman.diag(chains)
    ESS = effectiveSize(chains)
    MCSE = apply(beta, 2, sd) / sqrt(ESS)

    cat("  |**************************************************| 100%\n")

    cat('\nGelman-Rubin statistic [note: values close to 1.0 indicate convergence, and larger than 1.0 lack thereof]:\n')
    print(GB)
    cat('\nEffective Sample Size of the chains [note: values close to or larger than 10,000 are recommended]:\n')
    print(ESS)
    cat('\nMonte Carlo Standard Error [note: values are interpreted on the scale of the parameter]:\n')
    print(MCSE)
    cat('\n')

  }

  if(diagnostics == T){

    diagnost = list(GB, ESS, MCSE)
    return(list(beta.results, gamma.results, delta.results, diagnost))

  } else{

    return(list(beta.results, gamma.results, delta.results))

  }

}

