# SCEDbayes
Powerful multilevel Bayesian analysis, based on the Bayesian statistical model introduced by Swaminathan et al. (2014), to compute intervention effects in SCED studies

## Description
This function computes an autoregression at lag 1 time series analysis with Bayesian estimates. It can be used in single-case experimental designs to examine the effect of the introduction of a treatment B after a baseline A and the retreat of such treatment, also called ABAB reversal designs.

## Usage
ABABmodel(y, P, s, model = "level", plots = TRUE)

## Arguments
### y
outcome variable
### P 
phase identifier
### s 
session identifier
### model 
the model to be fitted. If set to "level" (the default), the model calculates intercepts only. Set to "trend", the model calculates intercepts and slopes.
### plots 
whether graphs are to be plotted. Defaults to TRUE.

## Return Values
### beta 
model coefficient estimates for A1, A1B1, B1A2, and A2B2 phase changes
### phases 
intercept and slope estimates for A1, B1, A2, and B2 phases
### delta 
effect size estimates for A1B1, B1A2, and A2B2 phase changes
