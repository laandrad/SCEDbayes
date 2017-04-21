## Plot time series
abab.its.plot.trend <- function(y, P, s, gamma) {

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)
  lPhase = c(0,nPhase, nPoints)
  lPhase = lPhase[-1] - lPhase[-length(lPhase)]
  lPhase = lPhase / 2

  plot(c(min(s)-1, max(s)+1), c(min(y)-1, max(y)+1), type='n', main='Interrupted Time Series Data Plot', frame.plot=F, ylab = 'Outcome', xlab = 'Session')
  points(y~s, type='b', pch=16, col='dark green')
  for(i in 1:3) abline(v=nPhase[i]+.5, lty = 2)
  nPhases = c(1, nPhase, length(s))
  phase.level = apply(gamma[,1:4], 2, mean)
  phase.slope = apply(gamma[,5:8], 2, mean)
  phases = c('A1', 'B1', 'A2', 'B2')
  mid.x = c(nPhases[1:4]+(nPhases[2:5]-nPhases[1:4])/2)
  text(mid.x, rep(max(y),4), phases, cex = 2)
  for(i in 1:4){
    x = c(nPhases[i]+.5, nPhases[i+1])
    y = c(phase.level[i], phase.level[i]+phase.slope[i]*lPhase[i])
    lines(x, y, lty = 2, col='red')
  }
  for(i in 1:4) text(nPhases[i], phase.level[i], bquote(paste(mu[.(phases[i])],' = ',.(round(phase.level[i],2)))), adj=c(-1,1))

}
