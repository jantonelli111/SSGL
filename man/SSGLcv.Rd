% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{SSGLcv}
\alias{SSGLcv}
\title{Find optimal lambda0 value using cross validation}
\usage{
SSGLcv(Y, X, lambda1, lambda0seq = seq(1, 100, by = 1), groups, a = 1,
  b = length(unique(groups)), nFolds = 5, M = 10, error = 0.001,
  forceGroups = c())
}
\arguments{
\item{Y}{The outcome to be analyzed}

\item{X}{An n by p matrix of covariates}

\item{lambda1}{Prior parameter for the slab component of the prior}

\item{lambda0seq}{Sequence of lambda0 values to consider}

\item{groups}{A vector of length p denoting which group each covariate is in}

\item{a}{First hyperparameter for the beta prior denoting the prior probability of being in the slab}

\item{b}{Second hyperparameter for the beta prior denoting the prior probability of being in the slab}

\item{nFolds}{The number of folds to run cross validation on}

\item{M}{Positive number less than p indicating how often to update theta and sigma. There is no
need to change this unless trying to optimize computation time}

\item{forceGroups}{A vector containing the indices of any groups you wish to automatically
include in the model and not penalize}
}
\value{
A list of values containing the lambda0 that minimizes the cross-validated error, the vector
        of cross validated errors for each lambda0 in the sequence, and the lambda0 sequence looked at
}
\description{
This function takes in an outcome and covariates, and uses cross validation
to find the best lambda0 value.
}
\examples{

## Here we generate 200 samples from 100 covariates
n = 200
G = 100
x = mvtnorm::rmvnorm(n, sigma=diag(G))


X = matrix(NA, nrow=n, ncol=G*2)
for (g in 1 : G) {
  X[,2*(g-1) + 1] = x[,g]
  X[,2*g] = x[,g]^2
}


Y = x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)

## Now find the best lambda0 using cross-validation
modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=.1, 
lambda0seq = seq(1,100, by=2),
groups = rep(1:G, each=2),
nFolds = 5)

modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=modSSGLcv$lambda0, 
               groups = rep(1:G, each=2))

modSSGL
}
