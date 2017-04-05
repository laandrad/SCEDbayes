## Plot time series
its.plot <- function(y, P, s, gamma) {

  nPoints = length(P)
  x = P[-1] - P[-nPoints]
  nPhase = which(!x==0)

  plot(y~s, main='Interrupted Time Series Data Plot', pch=16, type='b', col='dark green', frame.plot=F)
  for(i in 1:3) abline(v=nPhase[i]+.5, lty = 2)
  nPhases = c(1, nPhase, length(s))
  phase.level = apply(gamma[,1:4], 2, mean)
  phases = c('A1', 'B1', 'A2', 'B2')
  mid.x = c(nPhases[1:4]+(nPhases[2:5]-nPhases[1:4])/2)
  text(mid.x, rep(max(y),4), phases, cex = 2)
  for(i in 1:4) lines(c(nPhases[i]+.5, nPhases[i+1]+.5), rep(phase.level[i],2), lty = 2, col='red')
  for(i in 1:4) text(nPhases[i]+1, phase.level[i], bquote(paste(mu[.(phases[i])],' = ',.(phase.level[i]))), adj=c(-1,1))

}
