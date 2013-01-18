\name{Likelihood displacements}
\alias{recalculateLogLik}
\alias{obsDisp}
\alias{groupDisp}
\title{
Log-likelihood displacements for single observation and single grouping variable
}
\description{
Functions for log-likelihood displacements for each observation or each level of given factor
}
\usage{
recalculateLogLik(model, fixef = model@fixef, vcor = VarCorr(model)) 

groupDisp(formula, data, var) 

obsDisp(formula, data, inds=1:nrow(data)) 
}

\arguments{
  \item{model}{a mixed model of the class mer,}
  \item{fixef, vcor}{model parameters log-likelihood evaluation, if not provided then the estimates extracted from the 'model' parameter will be used}
  \item{formula}{a model formula that will be passes to the nlme function}
  \item{data}{a data frame}
  \item{var}{a name of grouping variable (factor) for which the group log-likelihood displacement will be performed}
  \item{inds}{indexes of observations for which observation log-likelihood displacement will be performed}
}
\details{
Likelihood displacement is defined as a difference of likelihoods calculated on full dataset for two models with different sets of parameters.  The first model is a model with ML estimates obtained for full dataset, while the second model is a model with ML estimates obtained on dataset without a selected observation or group of observations.

Likelihood displacements are used in model diagnostic.

Note that these functions reestimate coefficients in a set of model may be a time consuming.

The function recalculateLogLik() calculated a log-likelihood for model defined by the object model and model parameters defined in following function arguments.

The functions groupDisp() and obsDisp() calculates how the log-likelihood will decrees if selected groups or selected observations will not be used for parameter estimates. Note that log-likelihood is calculated on full dataset.

}
\author{
Przemyslaw Biecek
}
\examples{
data(eunomia)
require(lme4)
set.seed(1313)
eunomias <- eunomia[sample(1:2000,100),]
groupDisp(formula = BPRS.T2~ (1|CENTRE13), data=eunomias, var="CENTRE13")

obsDisp(formula = BPRS.T2~ (1|CENTRE13), data=eunomias, inds = 1:10)

obsDisp(formula = BPRS.T2~ (1|CENTRE13), data=eunomias)

}
