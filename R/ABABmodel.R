#' Bayesian Estimation of Effect Sizes in Single-Case Experimental Designs
#'
#' This function allows you to compute an interrupted time series (ITS) analysis with Bayesian estimates. It can be used in single-case experimental designs to examine the effect of the introduction of a treatment B after a baseline A and the retreat of such treatment thus the name of the model ABAB.
#' @param y outcome variable
#' @param P phase identifier
#' @param s session identifier
#' @param model the model to be fitted. It can be set to "level" (the default), "trend", or "AR1". "level" model calculates intercepts only, "trend" model calculates intercepts and slopes, and "AR1" calculates intercepts and an autocorrelation parameter at lag 1.
#' @param plots whether graphs are to be plotted. Defaults to TRUE.
#' @return delta effect size estimates for A1B1, B1A2, and A2B2 phase changes
#' @export
#' @examples
#' library(SCEDbayes)
#' dat = data(LAMBERT, package = "SCEDbayes")
#' dat = subset(dat, dat$STUDENT==1)
#' model = ABABmodel(dat$DATA.POINT, dat$PHASE, dat$SESSION, model = 'level', plots = TRUE)

ABABmodel = function(y, P, s, model = 'level', plots = TRUE) {

  ## load JAGS
  if(!require(rjags)){
    install.packages("rjags")
    library(rjags)
  }

  ## load model as specified in model argument
  if(model == 'trend'){
    ITS = model.trend

    } else if(model == 'AR1'){
    ITS = model.AR1

    } else{
    ITS = model.level
  }
  writeLines(ITS, con='model.txt')

  ## calculate bayesian coefficients using a Montecarlo Markov Chain
  beta = bayesian.estimates(y, P, s, model)
  beta = data.frame(beta)
  attach(beta)

  ## calculate phase intercept and slope at each point in the mcmc
  cat('Computing phase parameters...\n')

  if(model == 'trend'){
    gamma = phase.trends(beta)

    } else if(model == 'AR1'){
    gamma = phase.levels(beta)

    } else{
    gamma = phase.levels(beta)
  }
    cat("  |**************************************************| 100%\n")

  ## calculate effect size for each phase change according to model
  cat('Computing effect size...\n')

  if(model == 'trend'){
    delta = delta.trend(y, P, s, beta)

    } else if(model == 'AR1'){
    delta = delta.AR1(y, P, s, beta)

    } else{
    delta = delta.level(y, P, s, beta)
  }
  cat("  |**************************************************| 100%\n")

  ## Plotting results
  if(plots == TRUE){
    cat('Plotting results...\n')

    plot.titles = c('A1B1 Effect Size',
                    'B1A2 Effect Size',
                    'A2B2 Effect Size')

    parameter = c('delta A1B1',
                  'delta B1A2',
                  'delta A2B2')

    openGraph(width = 14, height = 7)

    if(model == 'trend'){
      layout(matrix(c(1:3,rep(4,3)), nrow = 2, byrow=T))
      sapply(1:3, function(i) posterior.plot(delta[,i], plot.titles[i], parameter[i]))
      its.plot.trend(y, P, s, gamma)

    } else{
      layout(matrix(c(1:3,rep(4,3)), nrow = 2, byrow=T))
      sapply(1:3, function(i) posterior.plot(delta[,i], plot.titles[i], parameter[i]))
      its.plot(y, P, s, gamma)

    }
    cat("  |**************************************************| 100%\n")

  } else{
    cat('Plots omitted...\n')

  }

  detach(beta)

  beta.results = sapply(beta, describe)
  gamma.results = sapply(gamma, describe)
  delta.results = sapply(delta, describe)

  print(beta.results)
  print(gamma.results)
  print(delta.results)

  return(list(beta, gamma, delta))

}

