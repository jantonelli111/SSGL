# SSGL
R package implementing the spike and slab group lasso

This is an R package to implement methods seen in "Spike-and-Slab Group Lassos for Grouped Regression and Sparse Generalized Additive Models" by Ray Bai, Gemma Moran, Joseph Antonelli, Yong Chen and Mary Boland, which can be found at the following link:

https://arxiv.org/pdf/1903.01979.pdf

To download the R package use the following in R:

```
library(devtools)
install_github(repo = "jantonelli111/SSGL")
library(SSGL)
```

# How to use the software

Here, we will simulate a simple example to show how the software works. First we will show how the function works for a chosen value of lambda0

```{r, eval=FALSE}
n = 200
G = 100
x = mvtnorm::rmvnorm(n, sigma=diag(G))


X = matrix(NA, nrow=n, ncol=G*2)

for (g in 1 : G) {
  X[,2*(g-1) + 1] = x[,g]
  X[,2*g] = x[,g]^2
}

Y = 200 + x[,1] + x[,2] + 0.6*x[,2]^2 + rnorm(n, sd=1)

## Now fit model for chosen lambda0 and lambda1 values
modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=10, 
groups = rep(1:G, each=2))

```

This will contain a list with estimates of the intercept, regression coefficients, residual variance, global sparsity parameter (theta), and the number of iterations it took to converge. Let's look at a couple of these now.

```{r, eval=FALSE}
## The first twenty regression coefficients
modSSGL$beta[1:20]
```
![Alt text](images/beta.png)

```{r, eval=FALSE}
## Estimate of the residual variance
modSSGL$sigmasq
```
![Alt text](images/sigmasq.png)

```{r, eval=FALSE}
## Estimate of theta
modSSGL$theta
```
![Alt text](images/theta.png)

So we see that our model has correctly identified the first two covariates (represented in this example by first 4 covariates) as nonzero, while zeroing out the remaining coefficients. This example relied on a well-chosen value of lambda0, which we won't know in general. To solve this, one can use the cross-validation function as below:


```{r, eval=FALSE}
modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=.1, lambda0seq = seq(4,20, by=2),
groups = rep(1:G, each=2), nFolds = 5)
```
![Alt text](images/CV1.png)

Notice that you might see a warning message about the algorithm diverging. This tends to happen only for small lambda0 values, and is not a problem unless you see this warning for each lambda0 value in your sequence. In our experience, lambda0 sequences should go from around 1 to 25. Now one can check what the chosen value of lambda0 is:

![Alt text](images/CV2.png)

And lastly run the model with the chosen lambda0 value:

```{r, eval=FALSE}
modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=modSSGLcv$lambda0, 
               groups = rep(1:G, each=2))
```

# Use in sparse generalized additive models (GAMs)

Now we will show how the SSGL procedure can be used to model sparse additive generalized models. We will use roughly the same functions as seen in the manuscript of the paper. We will model the effect of each covariate using 2 degree of freedom splines. First, we simulate G covariates, and then we create new design matrices using natural splines.

```{r, eval=FALSE}
library(splines)

TrueF = function(x) {
  return(5*sin(pi*x[,1]) + 
           2.5*(x[,3]^2 - 0.5) + 1*exp(x[,4]) + 3*x[,5])
}

n = 200
G = 100

## Generating training data
x = matrix(runif(G*n), nrow=n)
Y = TrueF(x) + rnorm(n, sd=1)

## degrees of freedom of the splines
mg = 2

X = matrix(NA, nrow=n, ncol=G*mg)

## create design matrix with 2 degree of freedom splines
## ensuring that the variables within each group are orthogonal
for (g in 1 : G) {
  splineTemp = ns(x[,g], df=mg)
  X[,mg*(g-1) + 1] = splineTemp[,1]
  for (m in 2 : mg) {
    tempY = splineTemp[,m]
    tempX = X[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]
    modX = lm(tempY ~ tempX)
    X[,mg*(g-1) + m] = modX$residuals
  }
}
```

And now we can build our SSGL model, first using cross-validation to find the appropriate tuning parameter value.

```{r, eval=FALSE}
lambda0seq = seq(3, 25, by=2)

## Cross validation
modSSGLcv = SSGLcv(Y=Y, X=X, lambda1=.1, 
                   lambda0seq = lambda0seq,
                   groups = rep(1:G, each=mg),
                   nFolds = 25)

## Final model
modSSGL = SSGL(Y=Y, X=X, lambda1=.1, lambda0=modSSGLcv$lambda0, 
               groups = rep(1:G, each=mg))
```

And finally, we can plot the effect of exposure 1 and exposure 3 on the outcome. 

```{r, eval=FALSE}
par(mfrow=c(1,2), pty='s')
## Plot the effect of exposure 1 on the outcome
ord = order(x[,1])
Curve = X[ord,1:2] %*% modSSGL$beta[1:2]
Truth = 5*sin(pi*x[ord,1])
plot(x[ord,1], Curve + mean(Truth) - mean(Curve), type='l', lwd=3,
     xlab="Covariate 1", ylab="f(X1)")
lines(x[ord,1], Truth, col=2, lwd=3)
legend("bottom", c("SSGL", "Truth"), lwd=3, lty=1, col=1:2)

## Plot the effect of exposure 3 on the outcome
ord = order(x[,3])
Curve = X[ord,5:6] %*% modSSGL$beta[5:6]
Truth = 2.5*(x[ord,3]^2 - 0.5)
plot(x[ord,3], Curve + mean(Truth) - mean(Curve), type='l', lwd=3,
     xlab="Covariate 3", ylab="f(X3)")
lines(x[ord,3], Truth, col=2, lwd=3)
legend("topleft", c("SSGL", "Truth"), lwd=3, lty=1, col=1:2)
```

![Alt text](images/Exposure1.png)
