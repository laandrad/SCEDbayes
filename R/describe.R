describe <- function(x) {
  estimate = c(mean(x), HDIofMCMC(x, .95))
  names(estimate) = c('estimate', 'lower.hdi', 'upper.hdi')
  return(estimate)
}
