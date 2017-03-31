# Example for Jags-Ymet-XmetSsubj-MrobustHier.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):

myData = read.csv( file="Lambert.csv" )
myData = cbind( myData , X3 = myData[,"X"] * myData[,"X2"] )
xName = "X" ; x2Name = "X2"; x3Name = "X3"; yName = "Y" ; sName="Subj"
fileNameRoot = "HierLinRegressData-Jags-" 

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XmetSsubj-MrobustHierExtendSCED.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , xName=xName ,x2Name = x2Name,x3Name = x3Name, yName=yName , sName=sName ,
                    numSavedSteps=20000 , thinSteps=15 , saveName=fileNameRoot )

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)

###################################################################################################################################################
# Start low level script
# Jags-Ymet-XmetSsubj-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
# This has the most updated version but won't work
source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , x2Name = "x2",x3Name = "x3", yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  x2 = data[,x2Name]
  x3 = data[,x3Name]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    x2 = x2,
    x3 = x3,
    y = y ,
    s = s ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
  Ntotal <- length(y)
  xm <- mean(x)
  x2m <- mean(x2)
  x3m <- mean(x3)
  ym <- mean(y)
  xsd <- sd(x)
  x2sd <- sd(x2)
  x3sd <- sd(x3)
  ysd <- sd(y)
  for ( i in 1:length(y) ) {
  zx[i] <- ( x[i] - xm ) / xsd
  zx2[i] <- ( x2[i] - x2m ) / x2sd
  zx3[i] <- ( x3[i] - x3m ) / x3sd
  zy[i] <- ( y[i] - ym ) / ysd
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx2[i] + zbeta3[s[i]] * zx3[i], 1/zsigma^2 , nu )
  }
  for ( j in 1:Nsubj ) {
  zbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )  
  zbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )
  zbeta2[j] ~ dnorm( zbeta2mu , 1/(zbeta2sigma)^2 )
  zbeta3[j] ~ dnorm( zbeta3mu , 1/(zbeta3sigma)^2 )
  
  
  }
  # Priors vague on standardized scale:
  zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta2mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta3mu ~ dnorm( 0 , 1/(10)^2 )
  zsigma ~ dnorm( 1.0E-3 , 1.0E+3 )
  zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta2sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta3sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  
  nu ~ dexp(1/30.0)
  # Transform to original scale:
  for ( j in 1:Nsubj ) {
  beta1[j] <- zbeta1[j] * ysd / xsd 
  beta2[j] <- zbeta2[j] * ysd / x2sd
  beta3[j] <- zbeta3[j] * ysd / x3sd
  beta0[j] <- zbeta0[j] * ysd  + ym - zbeta1[j] * xm * ysd / xsd + zbeta2[j] * x2m * ysd / x2sd + zbeta3[j] * x3m * ysd / x3sd
  }
  beta1mu <- zbeta1mu * ysd / xsd
  beta2mu <- zbeta2mu * ysd / x2sd
  beta3mu <- zbeta3mu * ysd / x3sd
  beta0mu <- zbeta0mu * ysd  + ym - zbeta1mu * xm * ysd / xsd + zbeta2mu * x2m * ysd / x2sd + zbeta3mu * x3m * ysd / x3sd
  sigma <- zsigma * ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1","beta2" , "beta3", "beta0mu" , "beta1mu" , "beta2mu", "beta3mu",
                  "zbeta0" , "zbeta1" , "zbeta2", "zbeta3",  "zbeta0mu" , "zbeta1mu" ,"zbeta2mu", "zbeta3mu",
                  "zsigma", "sigma", "nu" , "zbeta0sigma" , "zbeta1sigma","zbeta2sigma", "zbeta3sigma")
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}
