#' Estimate grouped spike and slab lasso regression model
#'
#' This function takes in an outcome and covariates, and estimates the 
#' posterior mode of the spike and slab group lasso penalty
#'
#' @param Y              The outcome to be analyzed
#' @param X              An n by p matrix of covariates
#' @param lambda1        Prior parameter for the slab component of the prior
#' @param lambda0        Prior parameter for the spike component of the prior
#' @param groups         A vector of length p denoting which group each covariate is in
#' @param a              First hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param b              Second hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param UpdateSigma    True or False indicating whether to estimate the residual variance
#' @param M              Positive number less than p indicating how often to update theta and sigma. There is no
#'                       need to change this unless trying to optimize computation time
#' @param betaStart      Starting values for beta vector. There is no need to change this unless for some specific reason
#' @param sigmaStart     Starting value for sigma. There is no need to change this unless for some specific reason
#' @param forceGroups    A vector containing the indices of any groups you wish to automatically
#'                       include in the model and not penalize
#'
#' @return A list of values containing the regression coefficients, the intercept, the estimate
#'         of the residual variance (this is simply the marginal variance of Y if UpdateSigma=FALSE),
#'         An estimate of theta the prior probability of entering into the slab, and the number of 
#'         iterations it took to converge 
#'
#' @export
#' @examples
#'
#' ## Here we generate 200 samples from 100 covariates
#' n = 200
#' G = 100
#' x = mvtnorm::rmvnorm(n, sigma=diag(G))
#' 
## Now create matrix that has linear and squared functions
## of the original covariates. 
#' 
#' X = matrix(NA, nrow=n, ncol=G*2)
#' 
#' for (g in 1 : G) {
#'   X[,2*(g-1) + 1] = x[,g]
#'   X[,2*g] = x[,g]^2
#' }
#' 
#' Y = 200 + x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)
#' 
#' ## Now fit model for chosen lambda0 and lambda1 values
#' modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=10, 
#' groups = rep(1:G, each=2))
#' 
#' modSSGL


SSGL = function(Y, X, lambda1, lambda0, groups,
                a = 1, b = length(unique(groups)),
                updateSigma = TRUE,
                M = 10, error = 0.001,
                betaStart = rep(0, dim(X)[2]),
                sigmasqStart = as.numeric(var(Y)),
                forceGroups = c()) {
  
  
  ## Number of groups and covariates overall
  G = length(unique(groups))
  p = length(groups)
  n = length(Y)
  
  ## Need to standardize them here
  means = apply(X, 2, mean)
  sds = apply(X, 2, sd)
  for (j in 1 : p) {
    X[,j] = (X[,j] - means[j]) / (sds[j]*sqrt(n-1))
  }
  
  ## Initialize values
  sigmasq = sigmasqStart
  beta = betaStart
  intercept = mean(Y - X %*% beta)
  Z = 1*((1:G) %in% unique(groups[which(beta != 0)]))
  theta = (a + sum(Z)) / (a + b + G)
  counter = 0
  diff = 10*error
  
  counter = 0
  
  while(diff > error & counter < 500) {
    ## Store an old beta so we can check for convergence at the end
    betaOld = beta
    
    ## First update intercept
    intercept = mean(Y - X %*% beta)
    
    ## Now update each group of parameters, beta_G
    for (g in 1 : G) {
      ## which parameters refer to this group
      active = which(groups == g)
      m = length(active)
      
      if (g %in% forceGroups) {
        yResid = Y - intercept - X[,-active] %*% as.matrix(beta[-active])
        beta[active] = solve(t(X[,active]) %*% X[,active]) %*% t(X[,active]) %*% yResid
      } else {
        
        ## Calculate delta for this size of a group
        gf = gFunc(beta = beta[active], lambda1 = lambda1,
                   lambda0 = lambda0, theta = theta, sigmasq = sigmasq)
        if (gf > 0) {
          delta =  sqrt(2*sigmasq*log(1/pStar(beta = rep(0,m), lambda1=lambda1,
                                              lambda0=lambda0, theta=theta))) +
            sigmasq*lambda1
        } else {
          delta = sigmasq*lambdaStar(beta = rep(0,m), lambda1=lambda1,
                                     lambda0=lambda0, theta=theta)
        }
        
        ## Calculate necessary quantities
        zg = t(X[,active]) %*% (Y - intercept - X[,-active] %*% as.matrix(beta[-active]))
        norm_zg = sum(zg^2)
        tempLambda = lambdaStar(beta = beta[active], lambda1 = lambda1,
                                lambda0 = lambda0, theta = theta)
        
        ## Update beta
        shrinkageTerm = (1 - sigmasq*lambdaStar(beta = beta[active], lambda1 = lambda1,
                                                lambda0 = lambda0, theta = theta)/norm_zg)
        shrinkageTerm = shrinkageTerm*(1*(shrinkageTerm > 0))
        beta[active] = shrinkageTerm*zg*(1*(norm_zg > delta))
      }
      
      
      ## Update Z
      Z[g] = 1*(any(beta[active] != 0))
      
      diff = sum((beta - betaOld)^2)
      
      if (g %% M == 0) {
        ## Update theta
        if (length(forceGroups) == 0) {
          theta = (a + sum(Z)) / (a + b + sum(Z)) 
        } else {
          theta = (a + sum(Z[-forceGroups])) / (a + b + sum(Z[-forceGroups])) 
        }
        
        ## Update sigmasq
        if (updateSigma) {
          sigmasq = sum((Y - X %*% beta - intercept)^2) / (n + 2)
          ## TODO do we need this? This avoids estimating sigma when it's too small
          ## and starts over the algorithm
          if (sigmasq < .00001*var(Y) | counter > 150) {
            sigmasq = var(Y)
            updateSigma = FALSE
            beta = rep(0, p)
            intercept = mean(Y)
            theta = 0.1
            diff=10*error
            Z = rep(0, G)
          }
        }
      }
    }
    
    ## Check to make sure algorithm doesn't explode for values of lambda0 too small
    tempSigSq = sum((Y - X %*% beta - intercept)^2) / (n + 2)
    if (tempSigSq < .00001*var(Y) | tempSigSq > 10000*var(Y)) {
      print(paste("Lambda0 = ", lambda0, " and algorithm is diverging. Increase lambda0"))
      sigmasq = var(Y)
      beta = rep(0, dim(X)[2])
      updateSigma = FALSE
      diff = 0
    }
    
    counter = counter + 1
    
  }
  
  ## need to re-standardize the variables
  betaSD = beta
  interceptSD = intercept
  for (j in 1 : p) {
    interceptSD = interceptSD - (beta[j]*means[j]) / (sds[j]*sqrt(n-1))
    betaSD[j] = beta[j] / (sds[j]*sqrt(n-1))
  }
  
  l = list(beta = betaSD, theta=theta, sigmasq=sigmasq, intercept=interceptSD, 
           nIter = counter)
  
  return(l)
}



#' Find optimal lambda0 value using cross validation
#'
#' This function takes in an outcome and covariates, and uses cross validation
#' to find the best lambda0 value.
#'
#' @param Y              The outcome to be analyzed
#' @param X              An n by p matrix of covariates
#' @param lambda1        Prior parameter for the slab component of the prior
#' @param lambda0seq     Sequence of lambda0 values to consider
#' @param nFolds         The number of folds to run cross validation on
#' @param groups         A vector of length p denoting which group each covariate is in
#' @param a              First hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param b              Second hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param UpdateSigma    True or False indicating whether to estimate the residual variance
#' @param M              Positive number less than p indicating how often to update theta and sigma. There is no
#'                       need to change this unless trying to optimize computation time
#' @param forceGroups    A vector containing the indices of any groups you wish to automatically
#'                       include in the model and not penalize
#'
#'
#' @return A list of values containing the lambda0 that minimizes the cross-validated error, the vector
#'         of cross validated errors for each lambda0 in the sequence, and the lambda0 sequence looked at
#'
#' @export
#' @examples
#'
#' ## Here we generate 200 samples from 100 covariates
#' n = 200
#' G = 100
#' x = mvtnorm::rmvnorm(n, sigma=diag(G))
#' 
## Now create matrix that has linear and squared functions
## of the original covariates.
#' 
#' X = matrix(NA, nrow=n, ncol=G*2)
#' for (g in 1 : G) {
#'   X[,2*(g-1) + 1] = x[,g]
#'   X[,2*g] = x[,g]^2
#' }
#' 
#' 
#' Y = x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)
#' 
#' ## Now find the best lambda0 using cross-validation
#' modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=.1, 
#' lambda0seq = seq(4,20, by=2),
#' groups = rep(1:G, each=2),
#' nFolds = 5)
#' 
#' modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=modSSGLcv$lambda0, 
#'                groups = rep(1:G, each=2))
#' 
#' modSSGL


SSGLcv = function(Y, X, lambda1, lambda0seq = seq(20, 50, by=1), 
                  groups, a = 1, b = length(unique(groups)),
                  nFolds = 5, updateSigma = TRUE,
                  M = 10, error = 0.001, forceGroups = c()) {
  
  NL = length(lambda0seq)
  n = length(Y)
  
  nVal = floor(n/nFolds)
  a11 = nFolds*(nVal + 1) - n
  holdouts = rep(1:nFolds, times=c(rep(nVal, a11), rep(nVal+1, nFolds - a11)))
  
  ## create grid storing cross validated error values
  gridError = matrix(NA, nrow=NL, ncol=nFolds)
  
  for (nf in 1 : nFolds)  {
    print(paste("Beginning fold # ", nf, sep=''))
    ## which subjects are excluded in the model estimation
    ws = which(holdouts == nf)
    
    ## Need to create new matrices which have norm 1
    Xtrain = X[-ws,]
    Ytrain = Y[-ws]
    ntrain = length(Ytrain)
    
    Xtest = X[ws,]
    Ytest = Y[ws]
    
    for (j in 1 : dim(X)[2]) {
      cmean = mean(Xtrain[,j])
      csd = sd(Xtrain[,j])
      
      Xtrain[,j] = (Xtrain[,j] - cmean) / (csd * sqrt(ntrain-1))
      Xtest[,j] = (Xtest[,j] - cmean) / (csd * sqrt(ntrain-1))
    }
    
    betaStart = rep(0, dim(Xtrain)[2])
    sigmasqStart = var(Ytrain)
    
    ## Loop through the different lambda0 values
    for (nl in 1 : NL) {
      lambda0 = lambda0seq[nl]
      modSSGL = SSGL(Y=Ytrain, X=Xtrain, lambda1=lambda1, lambda0=lambda0, 
                     groups = groups,
                     a = a, b = b,
                     updateSigma = updateSigma,
                     M = M, error = error,
                     betaStart = betaStart,
                     sigmasqStart = sigmasqStart,
                     forceGroups = forceGroups)
      
      betaStart = modSSGL$beta
      sigmasqStart = modSSGL$sigmasq
      
      yEst = modSSGL$intercept + Xtest %*% modSSGL$beta
      
      ## Keep track of squared error
      sqError = mean((yEst - Ytest)^2)
      gridError[nl,nf] = sqError
    }
  }
  
  ## average error across folds
  avgError = apply(gridError, 1, mean)
  
  ## minimum cross validated lambda
  lambda0 = lambda0seq[which(avgError == min(avgError))[1]]
  
  l = list(lambda0=lambda0,
           CVerror = avgError,
           lambda0seq = lambda0seq)
  
  return(l)
}
