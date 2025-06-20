require(skewt)
require(fitdistrplus)
source("Section2_functions.R")

alpha0 <- 0.1
alpha1 <- 0.3
beta1 <- 0.6
GARCHpars <- c(alpha0, alpha1, beta1)
set.seed(123)
sigma <- simgarch(dist = "norm", n=500000, pars = GARCHpars, copula = FALSE, outputsigma = TRUE)
summary(sigma)
# Calculate minimum sigma value
(k <- sqrt(alpha0/(1-beta1)))

# Find a parametric model for sigma
fitdata <- log(sigma - k)

#
dsnorm <- function (x, gamma = 1, mu = 0, sigma = 1, log = FALSE) 
{
  result <- rep(NaN, length(x))
  if ((sigma <= 0) | (gamma <= 0))
    return(result)
  x <- (x - mu)/sigma
  result[x < 0] <- dnorm(gamma * x[x < 0], log = TRUE)
  result[x >= 0] <- dnorm(x[x >= 0]/gamma, log = TRUE)
  result <- result + log(2/(gamma + 1/gamma)) - log(sigma)
  if (!(log))
    result <- exp(result)
  result
}

psnorm <- function (q, gamma = 1, mu = 0, sigma = 1) 
{
  result <- rep(NaN, length(q))
  if ((sigma <= 0) | (gamma <= 0))
    return(result)
  x <- (q - mu)/sigma
  result[x < 0] <- 2/(gamma^2 + 1) * pnorm(gamma * x[x < 0])
  result[x >= 0] <- 1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * 
    (pnorm(x[x >= 0]/gamma) - 1/2)
  result
}

qsnorm <- function (p, gamma = 1, mu = 0, sigma = 1) 
{
  result <- rep(NaN, length(p))
  if ((sigma <= 0) | (gamma <= 0))
    return(result)
  probzero <- 1/(gamma^2 + 1)
  result[p < probzero] <- 1/gamma * qnorm(((gamma^2 + 1) * p[p < 
                                                               probzero])/2)
  result[p >= probzero] <- gamma * qnorm((1 + 1/gamma^2)/2 * (p[p >= 
                                                                  probzero] - probzero) + 1/2)
  result * sigma + mu
}


# not sure about warnings in next code, but the fit seems fine
sigmamod <- fitdistrplus::fitdist(fitdata, "snorm", 
                                  start = list(gamma=0.9, mu=-1.5, sigma=0.75))
summary(sigmamod)
#check density
hist(fitdata, nclass = 60, prob=TRUE)
x <- seq(from = -5,to=5, length=500)
y <- dsnorm(x, gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3])
lines(x,y, col=2)
# check qqplot
qqplot(fitdata[1:10000], 
       qsnorm(ppoints(10000), gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3]),
       xlab="data",ylab="model")
abline(0,1)
abline(h=qsnorm(0.01, gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3]))
abline(v=qsnorm(0.01, gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3]))
abline(h=qsnorm(0.99, gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3]))
abline(v=qsnorm(0.99, gamma = coef(sigmamod)[1], mu = coef(sigmamod)[2], sigma = coef(sigmamod)[3]))

jointdens <- function(x,y){
  n1 <- length(x)
  n2 <- length(y)
  if ((n1 >1) & (n2==1)) y <- rep(y,n1)
  if ((n1 ==1) & (n2>1)) x <- rep(x,n2)
  n <- max(n1,n2)
  output <- rep(NA,n)
  innerfunc <- function(w,x,y){
    denom <- sqrt(alpha0 + alpha1*x^2 + beta1*w^2)
    dnorm(y/denom)*
      dnorm(x/w)*
      dsnorm(log(w-k), 
             gamma = coef(sigmamod)[1], 
             mu = coef(sigmamod)[2], 
             sigma = coef(sigmamod)[3])/
      (w*denom*(w-k))
  }
  for (i in 1:n)
  output[i] <- integrate(innerfunc, lower = k, upper = Inf, x=x[i], y=y[i])$value
  output
}

# Visualise joint density
ngrid <- 50
xvals <- seq(from = -3.5, to = 3.5, length=ngrid)
output <- outer(xvals, xvals, jointdens)
sum(output)*diff(xvals)[1]^2
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_GARCH_Gauss.pdf")
par(pty="s")
contour(xvals, xvals, output, nlevels=60)
abline(0,1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty = 3)
#dev.off()

# Find marginal density
margdens <- function(x){
  n <- length(x)
  output <- rep(NA,n)
  integrand <- function(w,x){
      dnorm(x/w)*
      dsnorm(log(w-k), 
             gamma = coef(sigmamod)[1], 
             mu = coef(sigmamod)[2], 
             sigma = coef(sigmamod)[3])/
      (w*(w-k))
  }
  for (i in 1:n)
    output[i] <- integrate(integrand, lower = k, upper = Inf, x=x[i])$value
  output
}
# Plot marginal density
curve(margdens, from = -4, to = 4)
integrate(margdens,-Inf,Inf)
# Find marginal cdf by numerical integration
margcdf <- function(x){
  n <- length(x)
  output <- rep(NA,n)
  for (i in 1:n)
    output[i] <- integrate(margdens, lower = -Inf, upper = x[i])$value
  output
}
curve(margcdf,from = -4, to = 4)
# Find marginal quantile function
margquant <- function(u){
  n <- length(u)
  output <- rep(NA,n)
  rootfunc <- function(x, u){
    margcdf(x) - u
  }
  for (i in 1:n)
    output[i] <- uniroot(rootfunc, interval = c(-10,10), u = u[i])$root
  output
}

# Copula density
ngrid <- 99
u <- seq(from = 0.01, to = 0.99, length=ngrid)
v <- u
xvals <- margquant(u)
fxvals <- margdens(xvals)
numerator <- outer(xvals,xvals,jointdens)
denominator <- outer(fxvals,fxvals)

#pdf("density_GARCH_Gauss.pdf")
par(pty="s")
contour(u, v, numerator/denominator, nlevels = 70, asp=1)
abline(0,1,lty=2)
abline(h = 0.5, lty=3)
abline(v = 0.5, lty=3)
#dev.off()

sum(numerator/denominator)*0.01^2
sum(abs(numerator/denominator -1))*0.01^2


# Remainder of code concerns conditional copula

conddens <- function(x, y, margin = 1){
  joint <- jointdens(x,y)
  switch(margin,
         joint/margdens(x),
         joint/margdens(y),
         stop("Invalid option"))
}

ngrid <- 50
xvals <- seq(from = -3.5, to = 3.5, length=ngrid)
# Conditional density of X_{t} given X_{t-1}
output <- outer(xvals,xvals,conddens,margin=1)
contour(xvals,xvals,output)
apply(output,1,sum)*diff(xvals)[1]
# Conditional density of X_{t-1} given X_{t}
output <- outer(xvals,xvals,conddens,margin=2)
contour(xvals,xvals,output)
apply(output,2,sum)*diff(xvals)[1]
# Compare condens
conddens(0.2,0.8, margin=1)
conddens(0.8,0.2, margin=2) # slightly different as expected

# Code to compute conditional cdf of X_{t} given X{t-1} or vice versa

cdf_integrand <- function(z,xory,margin=1){
  output <- rep(NA,length(z))
  innerfunc  <- function(w,x,y){
    denom <- sqrt(alpha0 + alpha1*x^2 + beta1*w^2)
    dnorm(y/denom)*
      dnorm(x/w)*
      dsnorm(log(w-k), 
             gamma = coef(sigmamod)[1], 
             mu = coef(sigmamod)[2], 
             sigma = coef(sigmamod)[3])/
      (w*denom*(w-k))
  }
  for (i in 1:length(output))
  output[i] <- switch(margin,
                      integrate(innerfunc, lower = k, upper = Inf, x=xory, y=z[i])$value,
                      integrate(innerfunc, lower = k, upper = Inf, x=z[i], y=xory)$value,
                      stop("Incorrect margin")
  )
  output
}



condcdf <- function(x, y, margin = 1){
  n1 <- length(x)
  n2 <- length(y)
  if ((n1 >1) & (n2==1)) y <- rep(y,n1)
  if ((n1 ==1) & (n2>1)) x <- rep(x,n2)
  n <- max(n1,n2)
  output <- rep(NA,n)
  for (i in 1:n)
  output[i] <- switch(margin,
    integrate(cdf_integrand, lower= -Inf, upper = y[i], xory = x[i], margin=1)$value,
    integrate(cdf_integrand, lower= -Inf, upper = x[i], xory = y[i], margin=2)$value,
    stop("Incorrect margin"))
  switch(margin,
         output/margdens(x),
         output/margdens(y),
         stop("Invalid margin"))
}

# Alternative method of calculation for X_t given X_{t-1}

condcdf1 <- function(x, y){
  n1 <- length(x)
  n2 <- length(y)
  if ((n1 >1) & (n2==1)) y <- rep(y,n1)
  if ((n1 ==1) & (n2>1)) x <- rep(x,n2)
  n <- max(n1,n2)
  output <- rep(NA,n)
  innerfunc <- function(w,x,y){
    denom <- sqrt(alpha0 + alpha1*x^2 + beta1*w^2)
    pnorm(y/denom)*
      dnorm(x/w)*
      dsnorm(log(w-k),
             gamma = coef(sigmamod)[1],
             mu = coef(sigmamod)[2],
             sigma = coef(sigmamod)[3])/
      (w*(w-k))
  }
  for (i in 1:n)
  output[i] <- integrate(innerfunc, lower = k, upper = Inf, x=x[i], y=y[i])$value/margdens(x[i])
 output
}


# Conditional cdf of X_{t} given X_{t-1}
output <- outer(xvals,xvals,condcdf,margin=1)
contour(xvals,xvals,output)
# Check on conditional cdf of X_{t} given X_{t-1}
condcdf(c(0.3, 0.7),c(0.9,0.5),margin=1)
condcdf1(c(0.3, 0.7),c(0.9,0.5))
# Conditional cdf of X_{t-1} given X_{t}
output <- outer(xvals,xvals,condcdf,margin=2)
contour(xvals,xvals,output)
condcdf(0.2, 0.9, margin=1)
condcdf(0.9, 0.2, margin=2) # slightly different as expected




# Conditional joint density of (X_{t-2}, X_{t}) given X_{t-1}
jointdens2 <- function(x, y, z){
  n1 <- length(x)
  n2 <- length(y)
  if ((n1 >1) & (n2==1)) y <- rep(y,n1)
  if ((n1 ==1) & (n2>1)) x <- rep(x,n2)
  n <- max(n1,n2)
  output <- rep(NA,n)
  innerfunc <- function(s, x, y, z){
    denom1 <- sqrt(alpha0 + alpha1*x^2 + beta1*s^2)
    denom2 <- sqrt(alpha0 + alpha1*z^2 + beta1*denom1^2)
    dnorm(y/denom2)*
      dnorm(z/denom1)*
      dnorm(x/s)*
      dsnorm(log(s-k), 
             gamma = coef(sigmamod)[1], 
             mu = coef(sigmamod)[2], 
             sigma = coef(sigmamod)[3])/
      (s*denom1*denom2*(s-k))
  }
  for (i in 1:n)
    output[i] <- integrate(innerfunc, lower = k, upper = Inf, x=x[i], y=y[i], z= z)$value
  #output/dst(z, df=coef(Xmod)[1], mu=coef(Xmod)[2], sigma =coef(Xmod)[3])
  output/margdens(z)
}
## Conditional density plot
ngrid <- 99
xvals <- seq(from = -3.5, to = 3.5, length=ngrid)
# Joint density of (X_{t-1}, X_{t+1}) given X_t
# middling value of X_t
output <- outer(xvals,xvals,jointdens2, z= 0)
contour(xvals,xvals,output)
sum(output)*diff(xvals)[1]^2
# more extreme value of X_t
xvals <- seq(from = -7, to = 7, length=ngrid)
output <- outer(xvals,xvals,jointdens2, z= 3)
contour(xvals,xvals,output,nlevels=30)
sum(output)*diff(xvals)[1]^2


# Final ingredient for conditional copula

# h-functions
hfunc <- function(u,v, margin=1){
  n1 <- length(u)
  n2 <- length(v)
  if ((n1 >1) & (n2==1)) v <- rep(v,n1)
  if ((n1 ==1) & (n2>1)) u <- rep(u,n2)
  x <- margquant(u)
  y <- margquant(v)
  output <- condcdf(x,
                    y,
                    margin)
  output <- pmax(pmin(output,1),0)
  output[(u==1) & (margin ==2)] <- 1
  output[(u==0) & (margin ==2)] <- 0
  output[(v==1) & (margin ==1)] <- 1
  output[(v==0) & (margin ==1)] <- 0
  output
}
# inverse h-functions
invhfunc <- function(u,v, margin=1){
  n1 <- length(u)
  n2 <- length(v)
  if ((n1 >1) & (n2==1)) v <- rep(v,n1)
  if ((n1 ==1) & (n2>1)) u <- rep(u,n2)
  n <- max(n1,n2)
  output <- rep(NA,n)
  rootfunc <- switch(margin,
                     function(x, u, v){hfunc(u, x, margin = 1) - v},
                     function(x, u, v){hfunc(x, v, margin = 2) - u},
                     stop("Invalid margin")) 
  for (i in 1:n){
    output[i] <- uniroot(rootfunc, interval = c(0.01,0.99),u = u[i], v = v[i])$root
  }
  output
}


####
# Checks on inverses
hfunc(c(0.3,0.9),c(0.1,0.8),margin=1)
invhfunc(c(0.3,0.9),hfunc(c(0.3,0.9),c(0.1,0.8),margin=1),margin=1)
plot(c(0.01,u,0.99), hfunc(0.3,c(0.01,u,0.99), margin=1))
hfunc(0.15,c(0.3,0.9),margin=2)
invhfunc(hfunc(0.15,c(0.3,0.9),margin=2),c(0.3,0.9),margin=2)
plot(c(0.01,u,0.99), hfunc(c(0.01,u,0.99),0.2, margin=2))
####


# Conditional copula plot
ngrid <- 99
u <- seq(from = 0.01, to = 0.99, length=ngrid)
v <- u
w <- 0.75
zval <- margquant(w)
xvals <- margquant(invhfunc(u,w,margin=2))
yvals <- margquant(invhfunc(w,v,margin=1))
numerator <- outer(xvals, yvals, jointdens2, z=zval)
fxvals <- conddens(xvals,zval,margin=2)
fyvals <- conddens(zval,yvals,margin=1)
denominator <- outer(fxvals,fyvals)
cop <- numerator/denominator
sum(cop)*diff(u)[1]^2
#pdf("density_GARCH_Gauss_cond_quartile.pdf")
par(pty="s")
contour(u, v, cop, nlevels = 70, asp=1)
abline(0,1,lty=2)
abline(h = 0.5, lty=3)
abline(v = 0.5, lty=3)
#dev.off()
sum(cop)*diff(u)[1]^2
sum(abs(cop -1))*diff(u)[1]^2


ngrid <- 99
u <- seq(from = 0.01, to = 0.99, length=ngrid)
v <- u
w <- 0.5
zval <- margquant(w)
xvals <- margquant(invhfunc(u,w,margin=2))
yvals <- margquant(invhfunc(w,v,margin=1))
numerator <- outer(xvals, yvals, jointdens2, z=zval)
fxvals <- conddens(xvals,zval,margin=2)
fyvals <- conddens(zval,yvals,margin=1)
denominator <- outer(fxvals,fyvals)
cop <- numerator/denominator
sum(cop)*diff(u)[1]^2
#pdf("density_GARCH_Gauss_cond_median.pdf")
par(pty="s")
contour(u, v, cop, nlevels = 70, asp=1)
abline(0,1,lty=2)
abline(h = 0.5, lty=3)
abline(v = 0.5, lty=3)
#dev.off()
sum(cop)*diff(u)[1]^2
sum(abs(cop -1))*diff(u)[1]^2

