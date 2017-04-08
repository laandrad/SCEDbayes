model.level = "
  model {
          # first observation
          y[1] ~ dnorm(mu[1], tau.e)
          mu[1] <- beta0 + beta1 * P1[1] + beta2 * P2[1] + beta3 * P3[1]

          # t subsequent observations
          for(t in 2:nPoints){
            y[t] ~ dnorm(mu[t], tau.delta)
            mu[t] <-  beta0 * (1-rho) +
                      beta1 * (P1[t]-rho*P1[t-1]) +
                      beta2 * (P2[t]-rho*P2[t-1]) +
                      beta3 * (P3[t]-rho*P3[t-1]) +
                      rho * y[t-1]
          }

          # Priors
          beta0 ~ dnorm(0, 0.001)
          beta1 ~ dnorm(0, 0.001)
          beta2 ~ dnorm(0, 0.001)
          beta3 ~ dnorm(0, 0.001)

          rho ~ dunif(-1,1)

          tau.e <- 1/pow(sigma.e,2)
          tau.delta <- 1/pow(sigma.delta,2)

          sigma.e <- sigma.delta/sqrt(1-pow(rho,2))
          sigma.delta ~ dunif(0,100)
        }
"

model.trend = "
  model {
          # first observation
          y[1] ~ dnorm(mu[1], tau.e)
          mu[1] <-  beta0         + beta4 * T1[1]         +  ## A1 intercept and slope
                    beta1 * P1[1] + beta5 * T2[1] * P1[1] +  ## A1B1 intercept and slope change
                    beta2 * P2[1] + beta6 * T3[1] * P2[1] +  ## B1A2 intercept and slope change
                    beta3 * P3[1] + beta7 * T4[1] * P3[1]    ## A2B2 intercept and slope change

          # t subsequent observations
          for(t in 2:nPoints){
            y[t] ~ dnorm(mu[t], tau.delta)
            mu[t] <-  beta0 * (1-rho)             + beta4 * (T1[t]       - rho*T1[t-1])         +  ## A1 intercept and slope
                      beta1 * (P1[t]-rho*P1[t-1]) + beta5 * (T2[t]*P1[t] - rho*T2[t-1]*P1[t-1]) +  ## A1B1 intercept and slope change
                      beta2 * (P2[t]-rho*P2[t-1]) + beta6 * (T3[t]*P2[t] - rho*T3[t-1]*P2[t-1]) +  ## B1A2 intercept and slope change
                      beta3 * (P3[t]-rho*P3[t-1]) + beta7 * (T4[t]*P3[t] - rho*T4[t-1]*P3[t-1]) +  ## A2B2 intercept and slope change
                      rho * y[t-1]
          }

          # Priors
          beta0 ~ dnorm(0, 0.001)
          beta1 ~ dnorm(0, 0.001)
          beta2 ~ dnorm(0, 0.001)
          beta3 ~ dnorm(0, 0.001)
          beta4 ~ dnorm(0, 0.001)
          beta5 ~ dnorm(0, 0.001)
          beta6 ~ dnorm(0, 0.001)
          beta7 ~ dnorm(0, 0.001)

          rho ~ dunif(-1,1)

          tau.e <- 1/pow(sigma.e,2)
          tau.delta <- 1/pow(sigma.delta,2)

          sigma.e <- sigma.delta/sqrt(1-pow(rho,2))
          sigma.delta ~ dunif(0,100)
        }
"
