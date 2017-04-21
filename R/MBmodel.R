#' Bayesian Estimation of Effect Sizes in Single-Case Experimental Designs
#'
#' This function allows you to compute an interrupted time series (ITS) analysis with Bayesian estimates. It can be used in single-case experimental designs to examine the effect of the introduction of a treatment B after a baseline A.
#' @param y outcome variable
#' @param P phase identifier
#' @param s session identifier
#' @param model the model to be fitted. If set to "level" (the default), the model calculates intercepts only. Set to "trend", the model calculates intercepts and slopes.
#' @param plots whether graphs are to be plotted. Defaults to TRUE.
#' @return delta effect size estimates for A1B1, B1A2, and A2B2 phase changes
#' @export
#' @examples
#' library(SCEDbayes)
#' dat = data(LAMBERT, package = "SCEDbayes")
#' dat = subset(dat, dat$STUDENT==1)
#' model = MBmodel(dat$DATA.POINT, dat$PHASE, dat$SESSION, model = 'level', plots = TRUE)

MBmodel = function(y, P, s, model = 'level', plots = TRUE, diagnostics = FALSE) {

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
  beta = mb.bayesian.estimates(y, P, s, model)
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

    if(diagnostics == T){
      openGraph(width = 7, height = 7)
      layout(1)
      gelman.plot(chains)
      ESS = effectiveSize(chains)
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
    ESS = effectiveSize(chains)
    cat('\nEffective Sample Size of the chains [note: values close to 10,000 are recommended]:\n')
    print(ESS)
  }

  return(list(beta.results, gamma.results, delta.results))

}

