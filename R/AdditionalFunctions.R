## Prior density
psi = function(beta, lambda) {
  m = length(beta)
  C = 2^(-m) * pi^(-(m-1)/2) / (gamma((m+1)/2))
  dens = C * lambda^m * exp(-lambda*sum(beta^2))
  
  return(dens)
}

## pStar function
pStar = function(beta, lambda1, lambda0, theta) {
  psi1 = psi(beta=beta, lambda=lambda1)
  psi0 = psi(beta=beta, lambda=lambda0)
  
  ## if a coefficient is really large then both these will 
  ## numerically be zero because R can't handle such small numbers
  if ((theta*psi1) == 0 & (1 - theta)*psi0 == 0) {
    p = 1
  } else {
    p = (theta*psi1) / (theta*psi1 + (1 - theta)*psi0)
  }
  
  return(p)
}

## Lambda star function
lambdaStar = function(beta, lambda1, lambda0, theta) {
  p = pStar(beta = beta, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)
  
  l = lambda1*p + lambda0*(1 - p)
  return(l)
}

## g function
gFunc = function(beta, lambda1, lambda0, theta, sigmasq) {
  l = lambdaStar(beta = beta, lambda1 = lambda1,
                 lambda0 = lambda0, theta = theta)
  p = pStar(beta = beta, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)
  
  g = (l - lambda1)^2 + (2/sigmasq)*log(p)
  return(g)
}

