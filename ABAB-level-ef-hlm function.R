# JAGS stands for Just Another Gibbs Sampler 
# JAGS builds MCMC sampler 
# JAGS takes a user's defined model and returns an MCMC sample of posterior distribution 
# To use the ABJags fuction, users first need to install JAGS and rjags 
# JAGS can be downloaded from http://mcmc-jags.sourceforge.net/

ABABJags = function(subj, session, y, P){
  if(!require(rjags)){
    install.packages("rjags")
    library(rjags)
  }
  
  model.s = "
  model {
    for(r in 1:nPoints){
      y[r] ~ dnorm(mu[r], tau.e[ subj[r] ])
      mu[r] <- beta0[ subj[r] ] + beta1[ subj[r] ] * P1[r] + beta2[ subj[r] ] * P2[r] + beta3[ subj[r] ] * P3[r]
    }

    # Prior specification for each subject
    for(s in 1:nSubj){
      beta0[s] ~ dnorm(beta0G, tau.beta0G)     
      beta1[s] ~ dnorm(beta1G, tau.beta1G)
      beta2[s] ~ dnorm(beta2G, tau.beta2G)     
      beta3[s] ~ dnorm(beta3G, tau.beta3G)

      tau.e[s] <- 1/pow(sigma.e[s],2)
      sigma.e[s] <- 1.0 / nu.sd[s]
      nu.sd[s] ~ dgamma(sG,rG)

    }

    # Hyper-priors on regression coefficients
    beta0G ~ dnorm(0, 0.001)     
    beta1G ~ dnorm(0, 0.001)
    beta2G ~ dnorm(0, 0.001)     
    beta3G ~ dnorm(0, 0.001)

    tau.beta0G ~ dgamma(.1, .1)
    tau.beta1G ~ dgamma(.1, .1)
    tau.beta2G ~ dgamma(.1, .1)
    tau.beta3G ~ dgamma(.1, .1)

    sG <- pow(m,2)/pow(d,2)
    rG <- m/pow(d,2)
    m ~ dgamma(.1,.1)
    d ~ dgamma(.1,.1)
  }
  "
  writeLines(model.s, con='model.txt')
  
  # create Phase dummy matrix
  max.n = max(session)
  m.p = matrix(NA, nrow=length(unique(subj)), ncol=max.n, byrow=T)
  for(j in 1:length(y)){
    s = subj[j]
    i = session[j]
    m.p[s,i] = P[j]
  }
  
  # translate matrix into vector form
  P1 = NULL
  P2 = NULL
  P3 = NULL
  for(s in 1:nrow(m.p)){
    phase = m.p[s,]
    phase = phase[!is.na(phase)]
    x = vector()
    for(i in 2:length(phase)) x = c(x, phase[i]-phase[i-1])
    nPhase = which(!x==0)
    P1 = c(P1, rep(0, nPhase[1]), rep(1, length(phase) - nPhase[1]))
    P2 = c(P2, rep(0, nPhase[2]), rep(1, length(phase) - nPhase[2]))
    P3 = c(P3, rep(0, nPhase[3]), rep(1, length(phase) - nPhase[3]))
  }
  
  # read number of points
  nPoints = length(y)
  t = 1:nPoints
  
  # read number of subjects
  nSubj = length(unique(subj))
  
  # put objects in list
  data = list(
    nPoints = nPoints,
    nSubj = nSubj,
    y = y,
    subj = subj,
    P1 = P1,
    P2 = P2,
    P3 = P3
  )
  
  # Initialize the chains with Maximum Likelihood estimators
  # Compute HLM with REML estimation method
  if(!require(lme4)){
    install.packages("lme4")
    library(lme4)
  }

  reg = lmer(y ~ P1 + P2 + P3 + (1|subj), REML = T, verbose = T)      #fit linear model
  reg = summary(reg)
  
  # Initial values for MCMC
  initsList = list(
    beta0G = reg$coefficients[1] ,    
    beta1G = reg$coefficients[2] ,       
    beta2G = reg$coefficients[3] ,    
    beta3G = reg$coefficients[4] 
  )

     # Specify Bayesian estimates to be extracted
  param = c('beta0', 'beta1', 'beta2', 'beta3', 'beta0G', 'beta1G', 'beta2G', 'beta3G')
  
  # Specify MCMC steps
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
  
  # summary(codaSamples)
  
  # Save the sampled chains as a matrix for further processing
  mcmcChain = as.matrix(codaSamples)
  
  # compute effect size
  # compute grand mean per phase (see Swaminathan et al (2014) p. 219 formula 15)
  phase = apply(cbind(P1,P2,P3),1,sum)
  
  grand.betas = mcmcChain[,c('beta0G','beta1G','beta2G','beta3G')]
  grand.means = apply(grand.betas, 1, cumsum)
  grand.means = t(grand.means)
  
  total.variance = function(y, grand.means, phase){
    grand.means = rep(grand.means, times=tapply(phase, phase, length))
    dev.sq.scores = apply(cbind(y,-grand.means), 1, function(y) sum(y)^2)
    sum.squares = tapply(dev.sq.scores, phase, sum)
    var.phase = sum.squares / tapply(phase, phase, length)
    var.tot = sum(var.phase)
    sigmaG = sqrt(var.tot) # sqrt of total variance = sigma e + sigma b
    return(sigmaG)
  }
  
  cat('Computing effect size...\n')
  sigmaG = apply(grand.means, 1, function(x) total.variance(y, x, phase))
  cat('  |**************************************************| 100% \n')
  
  n.beta1 = length(which(phase<2))
  n.beta3 = length(which(phase>1))
  
  ss.beta1 = mcmcChain[,'beta1G']*n.beta1
  ss.beta3 = mcmcChain[,'beta3G']*n.beta3
  delta = apply(cbind(ss.beta1,ss.beta3), 1, function(x) sum(x)/nPoints)
  delta = delta / sigmaG
  
  # Compute posterior central tendency and HDI
  posterior = data.frame(mcmcChain, delta)
  parameter = colnames(posterior)
  
  mu = apply(posterior, 2, mean)
  min.v = apply(posterior, 2, min)
  max.v = apply(posterior, 2, max)
  HDI = t(apply(posterior, 2, function(x) quantile(x, c(0.025,0.975))))
  holder = apply(posterior, 2, function(x) max(density(x)$y) / 6)
  
  # Graphing the results
  openGraph = function( width=7 , height=7 , ... ) {
    if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
      X11( width=width , height=height , type="cairo" , ... ) 
    } else { # Windows OS
      windows( width=width , height=height , ... )
    }
  }
  
  # plots
  myplot = function(i, label){
    # i = as.character(i)
    CD0 = HDI[i,]
    CD0 = ifelse(sum(CD0>0.2)==2 | sum(CD0 < -0.2)==2, '*** credibly not 0', '+ credibly 0')
    hdi.coords = xy.coords(HDI[i,], rep(holder[i],2))
    hist(posterior[,i], main=label, lwd=2, bty='n', probability = T, axes = F,
         sub = CD0, xlab='', ylab='', col='light green', border = 'white')
    axis(1, c(min.v, max.v), tck=0, labels = F)
    abline(v=c(-0.2,0.2), lty = 2, col='red')
    text(0, holder[i]*0.5, 'ROPE')
    lines(hdi.coords, lwd=3, col='dark green')
    text(HDI[i,],holder[i]*1.5,as.character(round(HDI[i,],2)))
    text(mu[i],holder[i]*5,bquote(paste(mu[.(as.character(i))],' = ',.(round(mu[i],2)))))
  }
  
  
  # plotting effect size
  openGraph(width = 14, height = 4)
  layout(matrix(1:3, ncol = 3, byrow=T))
  myplot('delta', 'Effect Size')
  myplot('beta1G', 'B1-A1 Level Difference')
  myplot('beta3G', 'B2-A2 Level Difference')
  
  
  # plotting time series
  y.p = split(y, subj)
  s.p = split(session, subj)
  Phase = apply(cbind(P1,P2,P3), 1, sum)
  p.p = split(Phase, subj)
  
  openGraph(width = 27, height = 13)
  ifelse(nSubj < 4, layout(matrix(1:3, ncol = 3, byrow=T)),
         ifelse(nSubj < 7, layout(matrix(1:6, ncol = 3, byrow=T)),
                ifelse(nSubj < 10, layout(matrix(1:9, ncol = 3, byrow=T)), 
                       layout(matrix(1:9, ncol = 3, byrow=T)))))
  n=min(nSubj,9)
  l = 2
  for(i in 1:n){
    y = y.p[[i]]
    s = s.p[[i]]
    p = p.p[[i]]
    plot(c(0,max(s)+l), c(0,max(y)+l), main=paste('Subject', i), type='n', frame.plot=F, xlab='session', ylab='outcome')
    lines(y~s, main=paste('Subject', i), type='b', pch=16, col='dark green')
    nPhase = tapply(p, factor(p), length)
    nPhase = cumsum(nPhase)
    for(j in 1:3) abline(v=nPhase[1:3]+.5, lty = 2)
    nPhases = c(1,s[nPhase])
    phase.level = mu[grep(x = names(mu), pattern = paste0('.', as.character(i),'.'), fixed = T)]
    phase.level = round(cumsum(phase.level),2)
    phases = c('A1', 'B1', 'A2', 'B2')
    mid.x = c(nPhases[1:4]+(nPhases[2:5]-nPhases[1:4])/2)
    text(mid.x, rep(max(y)+l,4), phases, cex = 1)
    for(j in 1:4) lines(c(nPhases[j]+.5, nPhases[j+1]+.5), rep(phase.level[j],2), lty = 2, col='red')
    # for(i in 1:4) text(nPhases[i]+1, phase.level[i], bquote(paste(mu[.(phases[i])],' = ',.(phase.level[i]))), adj=c(-1,1))
  } 
  
  # Outcome table
  results.table = data.frame(mu, HDI)
  colnames(results.table) = c('Estimate', 'Lower.HDI', 'Higher.HDI')
  
  group.level = results.table[grep('G', rownames(results.table)),]
  group.level = rbind(results.table['delta',], group.level)
  rownames(group.level) = c('Effect Size', rownames(group.level)[2:5])
  
  subj.level = results.table[-grep('G', rownames(results.table)),]
  subj.level = subj.level[-nrow(subj.level),]
  
  cat('\n', 'Table 1 \n',
      'Bayesian parameter estimates \n',
      'Subject Level \n')
  print(subj.level)
  
  cat('\n', 'Table 2 \n',
      'Bayesian parameter estimates \n',
      'Group Level \n')
  print(group.level)
  
}
