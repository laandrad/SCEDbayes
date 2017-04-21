## Plot time series
abab.its.plot <- function(y, P, s, gamma) {

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)

  plot(c(min(s)-1, max(s)+1), c(min(y)-1, max(y)+1), type='n', main='Interrupted Time Series Data Plot', frame.plot=F, ylab = 'Outcome', xlab = 'Session')
  points(y~s, type='b', pch=16, col='dark green')
  for(i in 1:3) abline(v=nPhase[i]+.5, lty = 2)
  nPhases = c(1, nPhase, length(s))
  phase.level = apply(gamma[,1:4], 2, mean)
  phases = c('A1', 'B1', 'A2', 'B2')
  mid.x = c(nPhases[1:4]+(nPhases[2:5]-nPhases[1:4])/2)
  text(mid.x, rep(max(y),4), phases, cex = 2)
  for(i in 1:4) lines(c(nPhases[i]+.5, nPhases[i+1]+.5), rep(phase.level[i],2), lty = 2, col='red')
  for(i in 1:4) text(nPhases[i]+1, phase.level[i], bquote(paste(mu[.(phases[i])],' = ',.(round(phase.level[i],2)))), adj=c(-1,1))

}
