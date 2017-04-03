# SCEDbayes
Powerful multilevel Bayesian analysis, based on the Bayesian statistical model introduced by Swaminathan et al. (2014), to compute intervention effects in SCED studies

## Description
This function allows you to compute an interrupted time series (ITS) analysis with Bayesian estimates. It can be used in single-case experimental designs to examine the effect of the introduction of a treatment B after a baseline A and the retreat of such treatment thus the name of the model ABAB.

## Usage
ABABmodel(y, P, s, model = "level", plots = TRUE)

## Arguments
y outcome variable
P phase identifier
s session identifier
model the model to be fitted. It can be set to "level" (the default), "trend", or "AR1". "level" model calculates intercepts only, "trend" model calculates intercepts and slopes, and "AR1" calculates intercepts and an autocorrelation parameter at lag 1.
plots whether graphs are to be plotted. Defaults to TRUE.

## Return Values
delta effect size estimates for A1B1, B1A2, and A2B2 phase changes
