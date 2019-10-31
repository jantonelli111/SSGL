#' Estimate grouped spike and slab lasso regression model along a path
#'
#' This function takes in an outcome and covariates, and estimates the 
#' posterior mode of the spike and slab group lasso penalty beginning at
#' a small value of lambda0 and continuing on a path until the chosen lambda0
#'
#' @param Y              The outcome to be analyzed
#' @param X              An n by p matrix of covariates
#' @param lambda1        Prior parameter for the slab component of the prior
#' @param lambda0        Prior parameter for the spike component of the prior
#' @param lambda0seq     Sequence of lambda0 grids to iterate through
#' @param groups         A vector of length p denoting which group each covariate is in
#' @param a              First hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param b              Second hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param M              Positive number less than p indicating how often to update theta and sigma. There is no
#'                       need to change this unless trying to optimize computation time
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
#' modSSGL = SSGLpath(Y=Y, X=X, lambda1=.1, lambda0=10, 
#' groups = rep(1:G, each=2))
#' 
#' modSSGL


SSGLpath = function(Y, X, lambda1, lambda0,
                    lambda0seq = seq(lambda1, lambda0, length=20),groups,
                    a = 1, b = length(unique(groups)),
                    M = 10, error = 0.001,
                    forceGroups = c()) {
  
  ## Final model
  betaStart = rep(0, dim(X)[2])
  updateSigma=FALSE
  printWarnings = FALSE
  for (nl in 1 : length(lambda0seq)) {
    lambda0 = lambda0seq[nl]
    
    if (nl == length(lambda0seq)) printWarnings = TRUE
    
    # starting values for lambda0 = lambda1
    if ( nl == 1) {
      modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                     groups = groups,
                     a = 1, b = G,
                     updateSigma = updateSigma,
                     M = 10, error = 0.001,
                     betaStart = betaStart,
                     theta = 0.5,
                     printWarnings = printWarnings
      )
      
    } else {
      modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                     groups = groups,
                     a = 1, b = G,
                     updateSigma = updateSigma,
                     M = 10, error = 0.001,
                     betaStart = betaStart,
                     sigmasqStart = sigmasqStart,
                     printWarnings = printWarnings)
      
    }
    
    
    if (modSSGL$nIter < 100 & modSSGL$converged) {
      updateSigma = TRUE
    }
    
    betaStart = modSSGL$betaStart
    sigmasqStart = modSSGL$sigmasqStart
    
  }
  
  l = list(beta = modSSGL$beta, theta=modSSGL$theta,
           sigmasq=modSSGL$sigmasq, intercept=modSSGL$intercept)
  
  return(l)
}



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
#' @param theta          Value of the sparsity parameter theta. There is no need to select a value for this parameter unless the 
#'                       sparsity is known a priori. If left blank, the function will update theta automatically 
#' @param printWarnings  Print warning messages about whether the optimization has converged
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
                sigmasqStart,
                theta,
                printWarnings = TRUE,
                forceGroups = c()) {
  
  ## Number of groups and covariates overall
  G = length(unique(groups))
  p = length(groups)
  n = length(Y)
  
  
  # get initial value for sigma
  df <- 3
  sigquant <- 0.9
  sigest <- sd(Y)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  
  if (missing(sigmasqStart)) {
    sigmasqStart <- sqrt(df * ncp / (df + 2))
  }
  
  Xstar = matrix(NA, dim(X)[1], dim(X)[2])
  for (j in 1 : dim(X)[2]) {
    Xstar[,j] = (X[,j] - mean(X[,j]))
  }
  
  ## store orthonormalized design matrix
  Xtilde = matrix(NA, dim(X)[1], dim(X)[2])
  
  ## store relevant matrices from SVD within each group
  Qmat = list()
  Dvec = list()
  
  ## create orthonormal matrix
  for (g in 1 : G) {
    active = which(groups == g)
    
    if (length(active) == 1) {
      Xtilde[,active] = sqrt(dim(Xstar)[1]) * (Xstar[,active] / sqrt(sum(Xstar[,active]^2)))
      Qmat[[g]] = NULL
      Dvec[[g]] = NULL
    } else {
      
      tempX = Xstar[,active]
      SVD = svd((1/n) * t(tempX) %*% tempX)
      Qmat[[g]] = SVD$u
      Dvec[[g]] = SVD$d
      
      Xtilde[,active] = Xstar[,active] %*% Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]]))
    }
  }
  
  converged=TRUE
  
  ## Initialize values
  sigmasq = sigmasqStart
  beta = betaStart
  intercept = mean(Y - Xtilde %*% beta)
  Z = 1*((1:G) %in% unique(groups[which(beta != 0)]))
  if (missing(theta)) {
    theta = (a + sum(Z))/(a + b + G)
  }
  counter = 0
  diff = 10*error
  
  counter = 0
  
  lambda0_base = lambda0
  
  while(diff > error & counter < 300) {
    
    ## Store an old beta so we can check for convergence at the end
    betaOld = beta
    
    ## Now update each group of parameters, beta_G
    for (g in 1 : G) {
      
      ## First update intercept
      active2 = which(beta != 0)
      if (length(active2) == 0) {
        intercept = mean(Y)
      } else if (length(active2) == 1) {
        intercept = mean(Y - Xtilde[,active2]*beta[active2])
      } else {
        intercept = mean(Y - Xtilde[,active2] %*% beta[active2]) 
      }
      
      ## which parameters refer to this group
      active = which(groups == g)
      m = length(active)
      lambda0 = sqrt(m) * lambda0_base
      
      if (g %in% forceGroups) {
        active2 = which(beta != 0 & !1:length(beta) %in% active)
        
        if (length(active2) == 0) {
          yResid = Y - intercept
        } else if (length(active2) == 1) {
          yResid = Y - intercept - Xtilde[,active2] * beta[active2]
        } else {
          yResid = Y - intercept - Xtilde[,active2] %*% as.matrix(beta[active2])
        }
        beta[active] = solve(t(Xtilde[,active]) %*% Xtilde[,active]) %*% t(Xtilde[,active]) %*% yResid
      } else {
        
        ## Calculate delta for this size of a group
        gf = gFunc(beta = rep(0, length(active)), lambda1 = lambda1,
                   lambda0 = lambda0, theta = theta, sigmasq = sigmasq, n=n)
        if (gf > 0) {
          delta =  sqrt(2*n*sigmasq*log(1/pStar(beta = rep(0,m), lambda1=lambda1,
                                                lambda0=lambda0, theta=theta))) +
            sigmasq*lambda1
        } else {
          delta = sigmasq*lambdaStar(beta = rep(0,m), lambda1=lambda1,
                                     lambda0=lambda0, theta=theta)
        }
        
        ## Calculate necessary quantities
        
        active2 = which(beta != 0 & !1:length(beta) %in% active)
        
        if (length(active2) == 0) {
          zg = t(Xtilde[,active]) %*% (Y - intercept)
        } else if (length(active2) == 1) {
          zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] * beta[active2])
        } else {
          zg = t(Xtilde[,active]) %*% (Y - intercept - Xtilde[,active2] %*% as.matrix(beta[active2]))
        }
        
        norm_zg = sqrt(sum(zg^2))
        tempLambda = lambdaStar(beta = beta[active], lambda1 = lambda1,
                                lambda0 = lambda0, theta = theta)
        
        ## Update beta
        shrinkageTerm = (1/n) * (1 - sigmasq*lambdaStar(beta = beta[active], lambda1 = lambda1,
                                                        lambda0 = lambda0, theta = theta)/norm_zg)
        shrinkageTerm = shrinkageTerm*(1*(shrinkageTerm > 0))
        beta[active] = as.numeric(shrinkageTerm)*zg*as.numeric((1*(norm_zg > delta)))
      }
      
      
      ## Update Z
      Z[g] = 1*(any(beta[active] != 0))
      
      diff = sqrt(sum((beta - betaOld)^2))
      
      if (g %% M == 0) {
        ## Update theta
        if (length(forceGroups) == 0) {
          theta = (a + sum(Z)) / (a + b + G) 
        } else {
          theta = (a + sum(Z[-forceGroups])) / (a + b + G - length(forceGroups)) 
        }
        
        ## Update sigmasq
        if (updateSigma) {
          active2 = which(beta != 0)
          if (length(active2) == 0) {
            sigmasq = sum((Y - intercept)^2) / (n + 2)
          } else if (length(active2) == 1) {
            sigmasq = sum((Y - Xtilde[,active2] * beta[active2] - intercept)^2) / (n + 2)
          } else {
            sigmasq = sum((Y - Xtilde[,active2] %*% as.matrix(beta[active2]) - intercept)^2) / (n + 2)
          }
        }
      }
    }
    
    ## Check to make sure algorithm doesn't explode for values of lambda0 too small
    active2 = which(beta != 0)
    if (length(active2) == 0) {
      tempSigSq = sum((Y - intercept)^2) / (n + 2)
    } else if (length(active2) == 1) {
      tempSigSq = sum((Y - Xtilde[,active2] * beta[active2] - intercept)^2) / (n + 2)
    } else {
      tempSigSq = sum((Y - Xtilde[,active2] %*% as.matrix(beta[active2]) - intercept)^2) / (n + 2)
    }
    
    if (updateSigma==FALSE & (tempSigSq < min_sigma2 | tempSigSq > 100*var(Y))) {
      if (printWarnings == TRUE) {
        print(paste("lambda0 = ", lambda0, ",", " Algorithm diverging. Increase lambda0 or lambda0seq", sep=""))
      }
      sigmasq = sigmasqStart
      betaStart = rep(0, dim(Xtilde)[2])
      converged = FALSE
      diff = 0
    }
    
    if (updateSigma==TRUE & (tempSigSq < min_sigma2 | tempSigSq > 100*var(Y))) {
      sigmasq = sigmasqStart
      beta = betaStart
      updateSigma = FALSE
      diff = 10*error
      counter = 0
    }
    
    if (sum(beta != 0) >= min(n, p)) {
      if (printWarnings == TRUE) {
        print("Beta is saturated. Increase lambda0 or lambda0seq")
      }
      sigmasq = sigmasqStart
      betaStart = rep(0, dim(Xtilde)[2])
      converged = F
      break
    }
    
    counter = counter + 1
    
  }
  
  ## need to re-standardize the variables
  betaSD = rep(NA, length(beta))
  for (g in 1 : G) {
    active = which(groups == g)
    if (length(active) == 1) {
      betaSD[active] = beta[active] * (sqrt(dim(Xstar)[1]) / sqrt(sum(Xstar[,active]^2)))
    } else {
      betaSD[active] = (Qmat[[g]] %*% diag(1/sqrt(Dvec[[g]])) %*% beta[active]) 
    }
  }
  
  ## Update new intercept
  interceptSD = mean(Y - X %*% betaSD)
  
  ## update the starting value for next iteration only if model converged
  if (converged == TRUE) {
    betaStart = beta
    if (updateSigma) {
      sigmasqStart = sigmasq
    }
  }
  
  
  ## estimate sigma regardless of convergence
  active2 = which(betaSD != 0)
  if (length(active2) == 0) {
    sigmasq = sum((Y - interceptSD)^2) / (n + 2)
  } else if (length(active2) == 1) {
    sigmasq = sum((Y - X[,active2] * betaSD[active2] - interceptSD)^2) / (n + 2)
  } else {
    sigmasq = sum((Y - X[,active2] %*% as.matrix(betaSD[active2]) - interceptSD)^2) / (n + 2)
  }
  
  l = list(beta = betaSD, betaStart = betaStart, theta=theta, sigmasqStart = sigmasqStart,
           sigmasq=sigmasq, intercept=interceptSD, nIter = counter,
           converged = converged)
  
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
#' lambda0seq = seq(1,100, by=2),
#' groups = rep(1:G, each=2),
#' nFolds = 5)
#' 
#' modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=modSSGLcv$lambda0, 
#'                groups = rep(1:G, each=2))
#' 
#' modSSGL

SSGLcv = function(Y, X, lambda1, lambda0seq = seq(1, 100, by=1), 
                  groups, a = 1, b = length(unique(groups)), nFolds = 5,
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
    
    # for (j in 1 : dim(X)[2]) {
    #   cmean = mean(Xtrain[,j])
    #   csd = sd(Xtrain[,j])
    # 
    #   Xtrain[,j] = (Xtrain[,j] - cmean)# / csd
    #   Xtest[,j] = (Xtest[,j] - cmean)# / csd
    # }
    
    sigmasqStart = var(Ytrain)
    betaStart = rep(0, dim(X)[2])
    updateSigma=FALSE
    printWarnings = FALSE
    ## Loop through the different lambda0 values
    for (nl in 1 : NL) {
      if (nl == NL) printWarnings = TRUE
      lambda0 = lambda0seq[nl]
      
      # starting values for lambda0 = lambda1
      if ( nl == 1) {
        modSSGL = SSGL(Y=Ytrain, X=Xtrain, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = 1,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       theta = 0.5)
        
      } else {
        modSSGL = SSGL(Y=Ytrain, X=Xtrain, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = 1,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       sigmasqStart = sigmasqStart)
        
      }
      
      
      if (modSSGL$nIter < 100 & modSSGL$converged) {
        updateSigma = TRUE
      }
      
      betaStart = modSSGL$betaStart
      sigmasqStart = modSSGL$sigmasqStart
      
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













#' Do sparse additive regression with splines using spike and slab group lasso
#'
#' This function takes in an outcome and covariates, and fits a semiparametric
#' regression model with additive functions so that the effect of each covariate is modeled
#' with a spline basis representation
#'
#' @param Y              The outcome to be analyzed
#' @param x              An n by p matrix of covariates
#' @param xnew           A new n2 by p matrix of covariates at which the outcome will be predicted
#' @param DF             Degrees of freedom of the splines for each covariate
#' @param lambda1        Prior parameter for the slab component of the prior
#' @param lambda0seq     Sequence of lambda0 values to consider
#' @param nFolds         The number of folds to run cross validation on
#' @param a              First hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param b              Second hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param UpdateSigma    True or False indicating whether to estimate the residual variance
#' @param error          The l2 norm difference that constitutes when the algorithm has converged
#' @param M              Positive number less than p indicating how often to update theta and sigma. There is no
#'                       need to change this unless trying to optimize computation time
#' @param forceGroups    A vector containing the indices of any groups you wish to automatically
#'                       include in the model and not penalize
#'
#'
#' @return A list of values containing the predictions at the observed locations, predictions
#'         at the new locations xnew, a vector indicating which covariates have nonzero coefficients, 
#'         an n by p matrix of predicted f(x) values for each covariate, an n2 by p matrix of 
#'         predicted f(xnew) values for each covariate, and an estimate of the residual
#'         variation in the model.
#' @export
#' @examples
#'
#' ## Here we generate 200 samples from 100 covariates
#' n = 200
#' n2 = 100
#' G = 100
#' x = mvtnorm::rmvnorm(n, sigma=diag(G))
#' xnew = mvtnorm::rmvnorm(n2, sigma=diag(G))
#' Ynew = 200 + xnew[,1] + xnew[,2] + 0.6*xnew[,2]^2 + rnorm(n2, sd=1)
#' 
#' modSSGLspr = SSGLspr(Y=Y, x=x, xnew = xnew, lambda1=.1, DF=2)

#' mean((modSSGLspr$predY - Y)^2)
#' mean((modSSGLspr$predYnew - Ynew)^2)
#' modSSGLspr$nonzero
#' 
#' plot(xnew[,2], modSSGLspr$fxnew[,2])
#' points(xnew[,2], xnew[,2] + 0.6*xnew[,2]^2 - 
#'          mean(xnew[,2] + 0.6*xnew[,2]^2), col=2)
#' legend("top", c("Estimated", "Truth"),col=1:2, pch=1)


SSGLspr = function(Y, x, xnew = NULL, DF=2, lambda1 = 1, 
                   lambda0seq = seq(1, 100, by=1), 
                   a = 1, b = dim(x)[2],
                   nFolds = 10, updateSigma = TRUE,
                   M = 10, error = 0.001, forceGroups = c()) {
  
  if (is.null(xnew) == FALSE) {
    
    mg = DF
    G = dim(x)[2]
    n = dim(x)[1]
    n2 = dim(xnew)[1]
    
    X = matrix(NA, nrow=n, ncol=G*mg)
    Xnew = matrix(NA, nrow=n2, ncol=G*mg)
    
    for (g in 1 : G) {
      splineTemp = splines::ns(x[,g], df=mg)
      splineTempNew = predict(splineTemp, xnew[,g])
      X[,mg*(g-1) + 1] = splineTemp[,1]
      Xnew[,mg*(g-1) + 1] = splineTempNew[,1]
      if (mg > 1) {
        for (m in 2 : mg) {
          tempY = splineTemp[,m]
          tempX = X[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]
          modX = lm(tempY ~ tempX)
          X[,mg*(g-1) + m] = modX$residuals
          Xnew[,mg*(g-1) + m] = splineTempNew[,m] - 
            cbind(rep(1,n2), Xnew[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]) %*% modX$coef
        }
      }
    }
    
    meansX = apply(X, 2, mean)
    sdX = apply(X, 2, sd)
    
    for (jj in 1 : dim(X)[2]) {
      X[,jj] = (X[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
      Xnew[,jj] = (Xnew[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
    }
    
    ## Cross validation
    modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=lambda1, 
                       lambda0seq = lambda0seq,
                       groups = rep(1:G, each=mg),
                       a = a, b = b, nFolds = nFolds,
                       M = M, error = error)
    
    ## Final model
    betaStart = rep(0, dim(X)[2])
    updateSigma=FALSE
    groups = rep(1:G, each=mg)
    printWarnings = FALSE
    for (nl in 1 : which.min(modSSGLcv$CVerror)) {
      lambda0 = lambda0seq[nl]
      
      if (nl == which.min(modSSGLcv$CVerror)) printWarnings=TRUE
      
      # starting values for lambda0 = lambda1
      if ( nl == 1) {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       theta = 0.5
        )
        
      } else {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       sigmasqStart = sigmasqStart)
        
      }
      
      
      if (modSSGL$nIter < 100 & modSSGL$converged) {
        updateSigma = TRUE
      }
      
      betaStart = modSSGL$betaStart
      sigmasqStart = modSSGL$sigmasqStart
      
    }
    
    predY = modSSGL$intercept + X %*% modSSGL$beta
    predYnew = modSSGL$intercept + Xnew %*% modSSGL$beta
    nonzero = 1*(modSSGL$beta[(1:G)*DF] != 0)
    sigmasq = modSSGL$sigmasq
    
    fx = matrix(NA, n, G)
    fxnew = matrix(NA, n2, G)
    
    groups = rep(1:G, each=DF)
    
    for (j in 1 : G) {
      if (DF > 1) {
        fx[,j] = as.vector(X[,which(groups == j)] %*% modSSGL$beta[which(groups == j)])
        fxnew[,j] = as.vector(Xnew[,which(groups == j)] %*% modSSGL$beta[which(groups == j)])
        fx[,j] = fx[,j] - mean(fx[,j])
        fxnew[,j] = fxnew[,j] - mean(fxnew[,j])
      } else {
        fx[,j] = as.vector(X[,which(groups == j)] * modSSGL$beta[which(groups == j)])
        fxnew[,j] = as.vector(Xnew[,which(groups == j)] * modSSGL$beta[which(groups == j)])
        fx[,j] = fx[,j] - mean(fx[,j])
        fxnew[,j] = fxnew[,j] - mean(fxnew[,j])
      }
    }
    
  } else {
    mg = DF
    G = dim(x)[2]
    n = dim(x)[1]

    X = matrix(NA, nrow=n, ncol=G*mg)

    for (g in 1 : G) {
      splineTemp = splines::ns(x[,g], df=mg)
      X[,mg*(g-1) + 1] = splineTemp[,1]
      if (mg > 1) {
        for (m in 2 : mg) {
          tempY = splineTemp[,m]
          tempX = X[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]
          modX = lm(tempY ~ tempX)
          X[,mg*(g-1) + m] = modX$residuals
        }
      }
    }
    
    meansX = apply(X, 2, mean)
    sdX = apply(X, 2, sd)
    
    for (jj in 1 : dim(X)[2]) {
      X[,jj] = (X[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
    }
    
    ## Cross validation
    modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=lambda1, 
                       lambda0seq = lambda0seq,
                       groups = rep(1:G, each=mg),
                       a = a, b = b, nFolds = nFolds,
                       M = M, error = error)
    
    ## Final model
    betaStart = rep(0, dim(X)[2])
    updateSigma=FALSE
    groups = rep(1:G, each=mg)
    printWarnings = FALSE
    for (nl in 1 : which.min(modSSGLcv$CVerror)) {
      lambda0 = lambda0seq[nl]
      
      if (nl == which.min(modSSGLcv$CVerror)) printWarnings = TRUE
      
      # starting values for lambda0 = lambda1
      if ( nl == 1) {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       theta = 0.5
        )
        
      } else {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       sigmasqStart = sigmasqStart)
        
      }
      
      
      if (modSSGL$nIter < 100 & modSSGL$converged) {
        updateSigma = TRUE
      }
      
      betaStart = modSSGL$betaStart
      sigmasqStart = modSSGL$sigmasqStart
      
    }
    
    predY = modSSGL$intercept + X %*% modSSGL$beta
    nonzero = 1*(modSSGL$beta[(1:G)*DF] != 0)
    sigmasq = modSSGL$sigmasq
    
    fx = matrix(NA, n, G)

    groups = rep(1:G, each=DF)
    
    for (j in 1 : G) {
      if (DF > 1) {
        fx[,j] = as.vector(X[,which(groups == j)] %*% modSSGL$beta[which(groups == j)])
        fx[,j] = fx[,j] - mean(fx[,j])
      } else {
        fx[,j] = as.vector(X[,which(groups == j)] * modSSGL$beta[which(groups == j)])
        fx[,j] = fx[,j] - mean(fx[,j])
      }
    }
    
    predYnew = NULL
    fxnew = NULL
  }
  
  
  l = list(predY = predY, predYnew = predYnew, nonzero = nonzero,
           fx = fx, fxnew = fxnew, sigmasq=sigmasq)
  
  return(l)
}






#' Fit nonlinear, 2-way interaction models with spike and slab group lasso
#'
#' This function takes in an outcome and covariates, and fits a semiparametric
#' regression model such that each covariate has a nonlinear main effect and
#' each pair of interaction terms has a nonlinear interaction term as well
#'
#' @param Y              The outcome to be analyzed
#' @param x              An n by p matrix of covariates
#' @param xnew           A new n2 by p matrix of covariates at which the outcome will be predicted
#' @param DFmain         Degrees of freedom of the splines for each main effect
#' @param DFint          Degrees of freedom of splines for each interaction pair. The total number of
#'                       terms in the interaction will be DFint*DFint
#' @param lambda1        Prior parameter for the slab component of the prior
#' @param lambda0seq     Sequence of lambda0 values to consider
#' @param nFolds         The number of folds to run cross validation on
#' @param a              First hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param b              Second hyperparameter for the beta prior denoting the prior probability of being in the slab
#' @param UpdateSigma    True or False indicating whether to estimate the residual variance
#' @param error          The l2 norm difference that constitutes when the algorithm has converged
#' @param M              Positive number less than p indicating how often to update theta and sigma. There is no
#'                       need to change this unless trying to optimize computation time
#' @param forceGroups    A vector containing the indices of any groups you wish to automatically
#'                       include in the model and not penalize
#'
#'
#' @return A list of values containing the predictions at the observed locations, predictions
#'         at the new locations xnew, a vector indicating which covariates have nonzero main effects,
#'         a matrix indicating which 2-way interactions are nonzero, and an estimate of the residual
#'         variation in the model.
#' @export
#' @examples
#'
#' ## Here we generate 200 samples from 100 covariates
#' n = 200
#' n2 = 100
#' G = 10
#' x = mvtnorm::rmvnorm(n, sigma=diag(G))
#' xnew = mvtnorm::rmvnorm(n2, sigma=diag(G))



SSGLint = function(Y, x, xnew = NULL, DFmain=2, DFint = 2,
                   lambda1 = 1, 
                   lambda0seq = seq(1, 100, by=1), 
                   a = 1, b = dim(x)[2],
                   nFolds = 10, updateSigma = TRUE,
                   M = 10, error = 0.001, forceGroups = c()) {
  
  if (is.null(xnew) == FALSE) {
    
    mgMain = DFmain
    mgInt = DFint
    
    G = dim(x)[2]
    n = dim(x)[1]
    n2 = dim(xnew)[1]
    
    X = matrix(NA, nrow=n, ncol=(G*mgMain) + choose(G,2)*(mgInt^2))
    Xnew = matrix(NA, nrow=n2, ncol=(G*mgMain) + choose(G,2)*(mgInt^2))
    
    for (g in 1 : G) {
      splineTemp = splines::ns(x[,g], df=mgMain)
      splineTempNew = predict(splineTemp, xnew[,g])
      X[,mgMain*(g-1) + 1] = splineTemp[,1]
      Xnew[,mgMain*(g-1) + 1] = splineTempNew[,1]
      if (mgMain > 1) {
        for (m in 2 : mgMain) {
          tempY = splineTemp[,m]
          tempX = X[,(mgMain*(g-1) + 1):(mgMain*(g-1) + m - 1)]
          modX = lm(tempY ~ tempX)
          X[,mgMain*(g-1) + m] = modX$residuals
          Xnew[,mgMain*(g-1) + m] = splineTempNew[,m] - 
            cbind(rep(1,n2), Xnew[,(mgMain*(g-1) + 1):(mgMain*(g-1) + m - 1)]) %*% modX$coef
        }
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        splineTemp1 = splines::ns(x[,g1], df=mgInt)
        splineTempNew1 = predict(splineTemp1, xnew[,g1])
        splineTemp2 = splines::ns(x[,g2], df=mgInt)
        splineTempNew2 = predict(splineTemp2, xnew[,g2])
        
        designX = matrix(NA, n, mgInt^2)
        designXnew = matrix(NA, n2, mgInt^2)
        
        counter2 = 1
        for (m1 in 1 : mgInt) {
          for (m2 in 1 : mgInt) {
            tempY = splineTemp1[,m1]*splineTemp2[,m2]
            tempMod = lm(tempY ~ splineTemp1 + splineTemp2)
            designX[,counter2] = tempMod$residuals
            tempYnew = splineTempNew1[,m1]*splineTempNew2[,m2]
            designXnew[,counter2] = tempYnew - cbind(rep(1,n2), splineTempNew1, splineTempNew2) %*%
              tempMod$coefficients
            counter2 = counter2 + 1
          }
        }
        
        X[,((G*mgMain) + (counter-1)*(mgInt^2) + 1) : ((G*mgMain) + counter*(mgInt^2))] = designX
        Xnew[,((G*mgMain) + (counter-1)*(mgInt^2) + 1) : ((G*mgMain) + counter*(mgInt^2))] = designXnew
        
        counter = counter + 1
      }
    } 
    
    meansX = apply(X, 2, mean)
    sdX = apply(X, 2, sd)
    
    for (jj in 1 : dim(X)[2]) {
      X[,jj] = (X[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
      Xnew[,jj] = (Xnew[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
    }
    
    groups = c(rep(1:G, each=mgMain), rep(G + 1:(choose(G,2)), each=mgInt^2))
    
    
    ## Cross validation
    modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=lambda1, 
                       lambda0seq = lambda0seq,
                       groups = groups,
                       a = a, b = b, nFolds = nFolds,
                       M = M, error = error)
    
    ## Final model
    betaStart = rep(0, dim(X)[2])
    updateSigma=FALSE
    printWarnings = FALSE
    for (nl in 1 : which.min(modSSGLcv$CVerror)) {
      lambda0 = lambda0seq[nl]
      
      if (nl == which.min(modSSGLcv$CVerror)) printWarnings = TRUE
      
      # starting values for lambda0 = lambda1
      if ( nl == 1) {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       theta = 0.5
        )
        
      } else {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       sigmasqStart = sigmasqStart)
        
      }
      
      
      if (modSSGL$nIter < 100 & modSSGL$converged) {
        updateSigma = TRUE
      }
      
      betaStart = modSSGL$betaStart
      sigmasqStart = modSSGL$sigmasqStart
      
    }
    
    predY = modSSGL$intercept + X %*% modSSGL$beta
    predYnew = modSSGL$intercept + Xnew %*% modSSGL$beta
    nonzero = 1*(modSSGL$beta[(1:G)*mgMain] != 0)
    sigmasq = modSSGL$sigmasq

    interactions = matrix(NA, G, G)

    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        interactions[g1,g2] = 1*(modSSGL$beta[((G*mgMain) + (counter-1)*(mgInt^2) + 1)] != 0)
        counter = counter + 1
      }
    }
    
  } else {
    mgMain = DFmain
    mgInt = DFint
    
    G = dim(x)[2]
    n = dim(x)[1]

    X = matrix(NA, nrow=n, ncol=(G*mgMain) + choose(G,2)*(mgInt^2))

    for (g in 1 : G) {
      splineTemp = splines::ns(x[,g], df=mgMain)
      X[,mgMain*(g-1) + 1] = splineTemp[,1]
      if (mgMain > 1) {
        for (m in 2 : mgMain) {
          tempY = splineTemp[,m]
          tempX = X[,(mgMain*(g-1) + 1):(mgMain*(g-1) + m - 1)]
          modX = lm(tempY ~ tempX)
          X[,mgMain*(g-1) + m] = modX$residuals
        }
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        splineTemp1 = splines::ns(x[,g1], df=mgInt)
        splineTemp2 = splines::ns(x[,g2], df=mgInt)

        designX = matrix(NA, n, mgInt^2)

        counter2 = 1
        for (m1 in 1 : mgInt) {
          for (m2 in 1 : mgInt) {
            tempY = splineTemp1[,m1]*splineTemp2[,m2]
            tempMod = lm(tempY ~ splineTemp1 + splineTemp2)
            designX[,counter2] = tempMod$residuals
            counter2 = counter2 + 1
          }
        }
        
        X[,((G*mgMain) + (counter-1)*(mgInt^2) + 1) : ((G*mgMain) + counter*(mgInt^2))] = designX

        counter = counter + 1
      }
    } 
    
    meansX = apply(X, 2, mean)
    sdX = apply(X, 2, sd)
    
    for (jj in 1 : dim(X)[2]) {
      X[,jj] = (X[,jj] - meansX[jj]) / (sdX[jj]* sqrt(n-1))
    }
    
    groups = c(rep(1:G, each=mgMain), rep(G + 1:(choose(G,2)), each=mgInt^2))
    
    
    ## Cross validation
    modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=lambda1, 
                       lambda0seq = lambda0seq,
                       groups = groups,
                       a = a, b = b, nFolds = nFolds,
                       M = M, error = error)
    
    ## Final model
    betaStart = rep(0, dim(X)[2])
    updateSigma=FALSE
    printWarnings = FALSE
    for (nl in 1 : which.min(modSSGLcv$CVerror)) {
      lambda0 = lambda0seq[nl]
      
      if (nl == which.min(modSSGLcv$CVerror)) printWarnings = TRUE
      
      # starting values for lambda0 = lambda1
      if ( nl == 1) {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       theta = 0.5
        )
        
      } else {
        modSSGL = SSGL(Y=Y, X=X, lambda1=lambda1, lambda0=lambda0, 
                       groups = groups,
                       a = 1, b = G,
                       updateSigma = updateSigma,
                       M = 10, error = 0.001,
                       betaStart = betaStart,
                       printWarnings = printWarnings,
                       sigmasqStart = sigmasqStart)
        
      }
      
      
      if (modSSGL$nIter < 100 & modSSGL$converged) {
        updateSigma = TRUE
      }
      
      betaStart = modSSGL$betaStart
      sigmasqStart = modSSGL$sigmasqStart
      
    }
    
    predY = modSSGL$intercept + X %*% modSSGL$beta
    predYnew = NULL
    nonzero = 1*(modSSGL$beta[(1:G)*mgMain] != 0)
    sigmasq = modSSGL$sigmasq
    
    interactions = matrix(NA, G, G)
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        interactions[g1,g2] = 1*(modSSGL$beta[((G*mgMain) + (counter-1)*(mgInt^2) + 1)] != 0)
        counter = counter + 1
      }
    }
  }
  
  
  l = list(predY = predY, predYnew = predYnew, nonzero = nonzero, 
           interactions = interactions, sigmasq = sigmasq)
  
  return(l)
}
