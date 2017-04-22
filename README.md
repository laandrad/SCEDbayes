# SCEDbayes
Powerful multilevel Bayesian analysis, based on the Bayesian statistical model introduced by Swaminathan et al. (2014), to compute intervention effects in SCED studies. The functions included in this package version use JAGS to sample the posterior space. The specification of the RJags algorithms are based on Kruschke (2010). The current version of the package includes two functions: MBmodel and ABABmodel. 

## Package Installation
The SCEDbayes package can be installed directly from GitHub into R using devtools package, as follows,

#### > install.packages("devtools")
#### > devtools::install_github("laandrad/SCEDbayes")
#### > library(SCEDbayes)

## ABABmodel
This function computes a lag-1 autoregression time series analysis with Bayesian estimates for ABAB reversal designs. It can be used with data from an ABAB reversal design to examine the effect of the introduction of a treatment B after a baseline A, followed up by a retreat and then a reintroduction of such treatment.

## MBmodel
This function computes a lag-1 autoregression time series analysis with Bayesian estimates for multiple-baseline designs. It can be used with data from a multiple-baseline design to examine the effect of the introduction of a treatment B after a baseline A.

## Data
The SCEDbayes package includes two data sets. The first set, called LAMBERT, is the Lambert et al. (2006) data set of an ABAB reversal design. The second set, called LAMBERT_MB, is a modified version of the first set which only includes the first two phases of the study per participant.

## Examples
### ABABmodel
> dat = subset(LAMBERT, LAMBERT$STUDENT==1)
> y = dat$DATA.POINT; P = dat$PHASE; s = dat$SESSION;
> model1 = ABABmodel(y, P, s, model = 'level') # Intercepts only model
> model2 = ABABmodel(y, P, s, model = 'trend') # Intercepts and slopes model

### MBmodel
> dat = subset(LAMBERT_MB, LAMBERT_MB$STUDENT==1)
> y = dat$DATA.POINT; P = dat$PHASE; s = dat$SESSION;
> model1 = MBmodel(y, P, s, model = 'level') # Intercepts only model
> model2 = MBmodel(y, P, s, model = 'trend') # Intercepts and slopes model

## References
Swaminathan, H., Rogers, H. J., & Horner, R. H. (2014). An effect size measure and Bayesian analysis of single-case designs. Journal of School Psychology, 52(2), 213-230. doi:10.1016/j.jsp.2013.12.002
Kruschke, J. (2010). Doing Bayesian Data Analysis: A Tutorial Introduction with R. Boston, MA.: Academic Press.
