mb.model.level = "
  model {
          # first observation
          y[1] ~ dt(mu[1], tau.e, nu)
          mu[1] <- beta0 + beta1 * P1[1]

          # t subsequent observations
          for(t in 2:nPoints){
            y[t] ~ dt(mu[t], tau.delta, nu)
            mu[t] <-  beta0 * (P1[t]-rho) +
                      beta1 * (P1[t]-rho*P1[t-1]) +
                      rho * y[t-1]
          }

          # Priors
          beta0 ~ dnorm(0, 0.001)
          beta1 ~ dnorm(0, 0.001)

          rho ~ dunif(-1,1)

          nu ~ dexp(1/30.0)

          tau.e <- 1/pow(sigma.e,2)
          tau.delta <- 1/pow(sigma.delta,2)

          sigma.e <- sigma.delta/sqrt(1-pow(rho,2))
          sigma.delta ~ dunif(0,100)
        }
"

mb.model.trend = "
  model {
          # first observation
          y[1] ~ dt(mu[1], tau.e, nu)
          mu[1] <-  beta0         + beta2 * T1[1]         +  ## A1 intercept and slope
                    beta1 * P1[1] + beta3 * T2[1] * P1[1]    ## A1B1 intercept and slope change

          # t subsequent observations
          for(t in 2:nPoints){
            y[t] ~ dt(mu[t], tau.delta, nu)
            mu[t] <-  beta0 * (1-rho)             + beta2 * (T1[t]       - rho*T1[t-1])         +  ## A1 intercept and slope
                      beta1 * (P1[t]-rho*P1[t-1]) + beta3 * (T2[t]*P1[t] - rho*T2[t-1]*P1[t-1]) +  ## A1B1 intercept and slope change
                      rho * y[t-1]
          }

          # Priors
          beta0 ~ dnorm(0, 0.001)
          beta1 ~ dnorm(0, 0.001)
          beta2 ~ dnorm(0, 0.001)
          beta3 ~ dnorm(0, 0.001)

          rho ~ dunif(-1,1)

          nu ~ dexp(1/30.0)

          tau.e <- 1/pow(sigma.e,2)
          tau.delta <- 1/pow(sigma.delta,2)

          sigma.e <- sigma.delta/sqrt(1-pow(rho,2))
          sigma.delta ~ dunif(0,100)
        }
"
