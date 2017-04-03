describe <- function(x) {
  estimate = c(mean(x), quantile(x, c(0.025, 0.975)))
  names(estimate) = c('estimate', 'lower.hdi', 'upper.hdi')
  return(estimate)
}
