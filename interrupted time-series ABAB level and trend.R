# JAGS stands for Just Another Gibbs Sampler 
# JAGS builds MCMC sampler 
# JAGS takes a user's defined model and returns an MCMC sample of posterior distribution 
# To use the ABJags fuction, users first need to install JAGS and rjags 
# JAGS can be downloaded from http://mcmc-jags.sourceforge.net/
# rjags can be installed by typing install.packages ("rjags") 

ABABJags = function(y,P){
  library(rjags)
  
  model.s = "
  model {
  for(i in 1:nPoints){
  y[i] ~ dnorm(mu[i], tau)    
  mu[i] <- beta0 + beta1 * T1[i] + beta2 * P1[i] + beta3 * T2[i] * P1[i] +
    beta4 * P2[i] + beta5 * T3[i] * P2[i] + beta6 * P3[i] + beta7 * T4[i] * P3[i]
  }
  beta0 ~ dnorm(0, 0.001)     
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  beta3 ~ dnorm(0, 0.001)
  beta4 ~ dnorm(0, 0.001)     
  beta5 ~ dnorm(0, 0.001)
  beta6 ~ dnorm(0, 0.001)
  beta7 ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
  }
  "
  writeLines(model.s, con='model.txt')
  
  # create dummy ABAB variables and Time according to Moeyaert et al. (2014)
  phase = P
  x = vector()
  for(i in 2:length(phase)) x = c(x, phase[i]-phase[i-1])
  nPhase = which(!x==0)
  P1 = c(rep(0, nPhase[1]), rep(1, length(phase)-nPhase[1]))
  P2 = c(rep(0, nPhase[2]), rep(1, length(phase)-nPhase[2]))
  P3 = c(rep(0, nPhase[3]), rep(1, length(phase) - nPhase[3]))
  T1 = 1:length(phase)
  T2 = T1 - (nPhase[1]+1)
  T3 = T1 - (nPhase[2]+1)
  T4 = T1 - (nPhase[3]+1)  
  y = y
  nPoints = length(y)
  # zy = (y - mean(y))/sd(y) # standardize data so that slope and intercept are not correlated
  
  data = list(
    nPoints = nPoints,
    y = y,
    P1 = P1,
    P2 = P2,
    P3 = P3,
    T1 = T1,
    T2 = T2,
    T3 = T3,
    T4 = T4
  )
  
  # Initialize the chains with Maximum Likelihood estimators
  reg = lm(y ~ T1 + P1 + T2*P1 + P2 + T3*P2 + P3 + T4*P3)      #fit linear model
  
  initsList = list(
    beta0 = reg$coefficients[1] ,    
    beta1 = reg$coefficients[2] ,       
    beta2 = reg$coefficients[3] ,       
    beta3 = reg$coefficients[4] ,       
    beta4 = reg$coefficients[5] ,    
    beta5 = reg$coefficients[6] ,       
    beta6 = reg$coefficients[7] ,       
    beta7 = reg$coefficients[8] ,       
    tau = length(y) / sum(reg$residuals^2)  
  )
  
  param = c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'beta5', 'beta6', 'beta7', 'tau')
  adaptSteps = 1000
  burnInSteps = 1000
  nChains = 3
  numSavedSteps = 20000
  thinSteps = 1
  nIter = ceiling((numSavedSteps * thinSteps)/nChains)
  
  # jagsModel = jags.model('model.txt', data=data, inits=initsList, n.chains = nChains, n.adapt = adaptSteps)
  jagsModel = jags.model('model.txt', data=data, n.chains = nChains, n.adapt = adaptSteps, inits=initsList)
  
  cat('Burning in the MCMC chain...\n')
  update(jagsModel, n.iter=burnInSteps)
  
  cat('Sampling final MCMC chain...\n')
  codaSamples = coda.samples(jagsModel, variable.names=param, n.iter=nIter, thin=thinSteps)
  
  summary(codaSamples)
  
  #Graph results
  mcmcChain = as.matrix(codaSamples)
  parameter = colnames(mcmcChain)
  
  return(mcmcChain)
  
  # windows()
  # layout(matrix(c(1,2,3,3), nrow = 2, byrow=T))
  # 
  # CD0 = NULL
  # 
  # for(i in 3:4){
  #   param = mcmcChain[,parameter[i]]
  #   param.mean = round(mean(param),3)
  #   holder = max(density(param)$y) / 6
  #   hdi = quantile(param, probs = c(0.025,0.975))
  #   CD0[i] = sum(hdi>0)
  #   hdi.coords = xy.coords(hdi, c(holder,holder))
  #   
  #   plot(density(param), main=bquote(paste('Posterior ',.(parameter[i]))), lwd=2, col='light green', bty='n')
  #   abline(v=0, lty = 2, col='red')
  #   lines(hdi.coords, lwd=3, col='blue')
  #   text(hdi[1],holder*2,as.character(round(hdi[1],2)))
  #   text(hdi[2],holder*2,as.character(round(hdi[2],2)))
  #   text(param.mean,holder*5,bquote(paste(mu,' = ',.(param.mean))))
  # }
  # 
  # CD0 = ifelse(CD0==1, 'not credibly different than 0', 'credibly different than 0')
  # 
  # # Phase A (Baseline)
  # A = subset(data.frame(t,y), P==0)
  # reg1 = lm(A$y ~ A$t)
  # 
  # # Phase B (Intervention)
  # B = subset(data.frame(t,y), P==1)
  # reg2 = lm(B$y ~ B$t)
  # 
  # holder = max(y) - max(y) / 20
  # subtitle = paste0('* Level is ', CD0[3], '; Slope is ', CD0[4])
  # 
  # # windows()
  # # layout(matrix(1), nrow = 1)
  # plot(y~t, main='Level and Trend Comparisons for Simulated Data', sub=subtitle)
  # abline(v=0, lty = 2)
  # lines(c(t[1],t[nA]), c(mean(A$y),mean(A$y)), lty = 2, col='red')
  # lines(c(t[nA+1],t[N]), c(mean(B$y),mean(B$y)), lty = 2, col='red')
  # lines(A$t, reg1$fitted.values, lty = 2, col='blue')
  # lines(B$t, reg2$fitted.values, lty = 2, col='blue')
  # text(c(t[3],t[length(t)-2]), c(holder, holder), c('Phase A', 'Phase B'))
  
}
