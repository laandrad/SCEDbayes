# JAGS stands for Just Another Gibbs Sampler 
# JAGS builds MCMC sampler 
# JAGS takes a user's defined model and returns an MCMC sample of posterior distribution 
# To use the ABJags fuction, users first need to install JAGS and rjags 
# JAGS can be downloaded from http://mcmc-jags.sourceforge.net/

ABABJags = function(y,P){
  if(!require(rjags)){
    install.packages("rjags")
    library(rjags)
  }
  
  model.s = "
  model {
  for(i in 1:nPoints){
  y[i] ~ dnorm(mu[i], tau)    
  mu[i] <- beta0 + beta1 * P1[i] + beta2 * P2[i] + beta3 * P3[i]
  }
  beta0 ~ dnorm(0, 0.001)     
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)     
  beta3 ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
  }
  "
  writeLines(model.s, con='model.txt')
  
  # create dummy ABAB variables and Time according to Moeyaert et al. (2014)
  phase = P
  x = vector()
  for(i in 2:length(phase)) x = c(x, phase[i]-phase[i-1])
  nPhase = which(!x==0)
  P1 = c(rep(0, nPhase[1]), rep(1, length(phase) - nPhase[1]))
  P2 = c(rep(0, nPhase[2]), rep(1, length(phase) - nPhase[2]))
  P3 = c(rep(0, nPhase[3]), rep(1, length(phase) - nPhase[3]))
  y = y
  nPoints = length(y)
  t = 1:nPoints

  data = list(
    nPoints = nPoints,
    y = y,
    P1 = P1,
    P2 = P2,
    P3 = P3
  )
  
  # Initialize the chains with Maximum Likelihood estimators
  reg = lm(y ~ P1 + P2 + P3)      #fit linear model
  
  initsList = list(
    beta0 = reg$coefficients[1] ,    
    beta1 = reg$coefficients[2] ,       
    beta2 = reg$coefficients[3] ,    
    beta3 = reg$coefficients[4] ,       
    tau = length(y) / sum(reg$residuals^2)  
  )
  
  param = c('beta0', 'beta1', 'beta2', 'beta3', 'tau') # Specify Bayesian estimates to be extracted
  adaptSteps = 1000
  burnInSteps = 1000
  nChains = 3
  numSavedSteps = 20000
  thinSteps = 1
  nIter = ceiling((numSavedSteps * thinSteps)/nChains)
  
  # Package everything to be sent to JAGS
  jagsModel = jags.model('model.txt', data=data, n.chains = nChains, n.adapt = adaptSteps, inits=initsList)
  
  # Start the chains and get rid of first n steps
  cat('Burning in the MCMC chain...\n')
  update(jagsModel, n.iter=burnInSteps)
  
  # Get the tunned chains
  cat('Sampling final MCMC chain...\n')
  codaSamples = coda.samples(jagsModel, variable.names=param, n.iter=nIter, thin=thinSteps)
  
  summary(codaSamples)
  
  # Save the sampled chains as a matrix for further processing
  mcmcChain = as.matrix(codaSamples)
  parameter = colnames(mcmcChain)

  # Graphing the results
  openGraph = function( width=7 , height=7 , ... ) {
    if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
      X11( width=width , height=height , type="cairo" , ... ) 
    } else { # Windows OS
      windows( width=width , height=height , ... )
    }
  }

  openGraph(width = 14, height = 7)
  layout(matrix(c(1:4,rep(5,4)), nrow = 2, byrow=T))

  plot.titles = c('A1 Level', 'B1-A1 Level Difference',
                  'A2-B1 Level Difference', 'B2-A2 Level Difference')
  
  mu = apply(mcmcChain, 2, mean)
  min.v = apply(mcmcChain, 2, min)
  max.v = apply(mcmcChain, 2, max)
  HDI = t(apply(mcmcChain, 2, function(x) quantile(x, probs = c(0.025,0.975))))
  holder = apply(mcmcChain, 2, function(x) max(density(x)$y) / 6)
  
  for(i in 1:4){
    CD0 = HDI[i,]>0
    CD0 = ifelse(CD0[1]==CD0[2], '*** credibly not 0', '+ credibly 0')
    hdi.coords = xy.coords(HDI[i,], rep(holder[i],2))
    hist(mcmcChain[,i], main=plot.titles[i], lwd=2, bty='n', probability = T, axes = F,
         sub = CD0, xlab='', ylab='', col='light green', border = 'white')
    axis(1, c(min.v, max.v), tck=0, labels = F)
    abline(v=0, lty = 2, col='red')
    text(0, holder[i]*1.5, '0')
    lines(hdi.coords, lwd=3, col='dark green')
    text(HDI[i,1],holder[i]*1.5,as.character(round(HDI[i,1],2)))
    text(HDI[i,2],holder[i]*1.5,as.character(round(HDI[i,2],2)))
    text(mu[i],holder[i]*5,bquote(paste(mu[.(parameter[i])],' = ',.(mu[i]))))
  }
  
  plot(y~t, main='Level Comparisons', pch=16, type='b', col='dark green', frame.plot=F)
  for(i in 1:3) abline(v=nPhase[i]+.5, lty = 2)
  nPhases = c(1, nPhase, length(t))
  phase.level = round(cumsum(mu),2)
  phases = c('A1', 'B1', 'A2', 'B2')
  mid.x = c(nPhases[1:4]+(nPhases[2:5]-nPhases[1:4])/2)
  text(mid.x, rep(max(y),4), phases, cex = 2)
  for(i in 1:4) lines(c(nPhases[i]+.5, nPhases[i+1]+.5), rep(phase.level[i],2), lty = 2, col='red')
  for(i in 1:4) text(nPhases[i]+1, phase.level[i], bquote(paste(mu[.(phases[i])],' = ',.(phase.level[i]))), adj=c(-1,1))
  
  results.table = data.frame(parameter[1:length(mu)], mu, HDI)
  colnames(results.table) = c('Parameter', 'Estimate', 'Lower.HDI', 'Higher.HDI')

  return(results.table)
  
}
