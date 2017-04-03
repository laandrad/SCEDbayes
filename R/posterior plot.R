## posterior distribution plot
posterior.plot <- function(x, title, parameter) {

  mu = mean(x)
  min.v = min(x)
  max.v = max(x)
  HDI = quantile(x, c(0.025,0.975))
  holder = max(density(x)$y) / 6
  CD0 = ifelse(sum(HDI>0.2)==2 | sum(HDI < -0.2)==2, '* credibly not 0', '+ credibly 0')
  hdi.coords = xy.coords(HDI, rep(holder,2))

  hist(x, main=title, lwd=2, bty='n', probability = T, axes = F,
       sub = CD0, xlab='', ylab='', col='light green', border = 'white')
  axis(1, c(min.v, max.v), tck=0, labels = F)
  abline(v=c(-0.2,0.2), lty = 2, col='red')
  text(0, holder*0.5, 'ROPE')
  lines(hdi.coords, lwd=3, col='dark green')
  text(HDI,holder*1.5,as.character(round(HDI,2)))
  text(mu,holder*5,bquote(paste(mu[.(parameter)],' = ',.(round(mu,2)))))

}
