library(tscopula)
library(rvinecopulib)
source("Dvine_functions.R")
# Data from LMSM18 paper:
load("LMSM18data.RData")

# These are the log-returns on the USD-AUD exchange rate:
ts.plot(Xt)
# These are the pseudo-observations from the copula process obtained in LMSM18
# by a kernel density estimate:
ts.plot(Ut)
# Should be approx uniform
hist(Ut)

# Copula T1-T5
# Fitted with tscopula package
(T1fit <- fit(dvinecopula(family = "t", pars = list(c(0, 10))), Ut))
(T2fit <- fit(dvinecopula(family = "t", pars = list(c(0, 7), c(0, 10))), Ut))
# Exact results for T5 are sensitive to starting values; slowish fit
(T5fit <- fit(dvinecopula(family = "t", 
                       pars = list(c(0, 6),c(0, 7), c(0, 8), c(0, 9), c(0, 10))), Ut, 
           control = list(maxit = 5000)))
AIC(T1fit, T2fit, T5fit)
BIC(T1fit, T2fit, T5fit)

# Copula A1-A5
uvals <- seq(from = 0.01, to = 0.99, length=50)
z <- outer(uvals, uvals, dtmix, theta = c(0.01, 5, 0.01, 5, 0.5))
contour(z)
z2 <- outer(uvals, uvals, htmix, theta = c(0.01, 5, 0.01, 5, 0.5), cond_var = 1)
contour(z2)
#
(A1fit <- fit_tmix(Ut, startpars = c(0.5, 6, 0.5, 4, 0.5)))
(A2fit <- fit_tmix(Ut, startpars = rep(c(0.01, 7, 0.01, 4, 0.8), 2), order = 2, maxit =2000))
# very slow fit:
thetaA0 <- c(0.0095, 26.3887, 0.0955, 2.5903, 0.7313, 
            0.0954, 30.2724, 0.0575, 6.2465, 0.1843, 
            0.0612, 37.8164, 0.0372, 6.3411, 0.0977, 
            0.0846, 10.0860, 0.7714, 4.4363, 0.9719, 
            0.2057, 27.8140, 0.3705, 25.0548, 0.6650) #the posterior values from MATLAB
system.time(A5fit <- fit_tmix(Ut, startpars = thetaA0, order = 5, maxit = 20000))
AIC(A1fit, A2fit, A5fit)
BIC(A1fit, A2fit, A5fit)

# Copula B1-B5
uvals <- seq(from = 0.01, to = 0.99, length=50)
z <- outer(uvals, uvals, dCGmix, theta = c(0.4, 0.7, 1.5, 1.2, 0.5))
contour(z)
z2 <- outer(uvals, uvals, hCGmix, theta = c(0.4, 0.7, 1.5, 1.2, 0.5), cond_var = 2)
contour(z2)
#
(B1fit <- fit_CGmix(Ut, startpars = c(0.5, 0.5, 2, 2, 0.5)))
(B2fit <- fit_CGmix(Ut, startpars = rep(c(0.3, 0.75, 1.2, 1.2, 0.5), 2), order = 2))
# very slow fit:
thetaB0 <- c(0.5305, 0.4083, 1.1857, 1.1492, 0.5113, 
            0.6774, 0.4070, 1.1562, 1.4388, 0.6571, 
            0.5801, 0.0696, 1.0897, 1.4539, 0.7347, 
            0.3893, 0.1708, 1.5307, 1.0585, 0.1732, 
            0.3015, 0.5624, 2.1533, 1.0277, 0.0515) #the posterior value from MATLAB
system.time(B5fit <- fit_CGmix(Ut, startpars = thetaB0, order = 5, maxit = 20000))
AIC(B1fit, B2fit, B5fit)
BIC(B1fit, B2fit, B5fit)

# Using inverse-v-transformed copulas

# AR(1)-type models
startpars <- list(ar = 0.3)
joeAR1 <- fit(dvinecopulavt(family = "joe", pars = startpars), Ut)
claytonAR1 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars), Ut)
astAR1 <- fit(dvinecopulavt(family = "ast", pars = startpars), Ut)

# AR(5)-type models
startpars <- list(ar = c(0.05, 0.05, 0.05, 0.05, 0.05))
# Joe and Clayton need increased iterations
joeAR5 <- fit(dvinecopulavt(family = "joe", pars = startpars), Ut)
claytonAR5 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars), Ut,
                  control= list(maxit = 10000))
astAR5 <- fit(dvinecopulavt(family = "ast", pars = startpars), Ut)

# ARMA(1,1)-type models
startpars <- list(ar = 0.95, ma = -0.85)
joeARMA11 <- fit(dvinecopulavt(family = "joe", pars = startpars, 
                               maxlag = 40), Ut)
claytonARMA11 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars,
                                   maxlag = 40), Ut)
astARMA11 <- fit(dvinecopulavt(family = "ast", pars = startpars,
                               maxlag = 40), Ut,
                 tsoptions= list(hessian = TRUE))

# Results
allmods <- list("T1" = T1fit, "A1" = A1fit, "B1" = B1fit, 
                "Joe-AR(1)" = joeAR1, "Clayton180-AR(1)" = claytonAR1, "ast-AR(1)" = astAR1,
                "T5" = T5fit, "A5" = A5fit, "B5" = B5fit,
                "Joe-AR(5)" = joeAR5, "Clayton180-AR(5)" = claytonAR5, "ast-AR(5)" = astAR5,
                "Joe-ARMA(1,1)" = joeARMA11, "Clayton180-ARMA(1,1)" = claytonARMA11, 
                "ast-ARMA(1,1)" = astARMA11)
# Convergence check
sapply(allmods, function(mod){
  if (is(mod, "tscopulafit"))
    output <- mod@fit$convergence
  else
    output <- mod$conv
  output})    
#
allAIC <- sapply(allmods, AIC)
allBIC <- sapply(allmods, BIC)
#
npars <- sapply(allmods, function(mod){
  if (is(mod, "tscopulafit"))
    output <- length(mod@fit$par)
  else
    output <- length(mod$par.ests)
  output})  
logLik <- round(-(allAIC - 2*npars)/2, 2)
#
order <- sapply(allmods, function(mod){
  if (is(mod, "tscopulafit")){
    if (is(mod@tscopula, "dvinecopulavt"))
      output <- mod@fit$EML
    else
      output <- length(mod@fit$par)/2
  }
  else
    output <- mod$npar/5
  output}) 
#
allresults <- data.frame("D-vine order" = order, "No pars" = npars, Loglik = round(logLik,1), AIC = round(allAIC,1), 
                         BIC = round(allBIC,1))
allresults
xtable::xtable(allresults, display = c("s","d", "d", "f", "f", "f"), digits =1)


astARMA11
(nu <- sapply(tscopula:::mklist_dvinevt(astARMA11@tscopula,2),
             function(cop){cop$parameters}))
kendall(astARMA11, 2)



# VaR prediction
alpha <- c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99)
n <- length(Ut)

# Methods to do: AR(5)-Joe, AR(5)-ast, ARMA(1,1)-clayton, ARMA(1,1)-ast (just best models)
# Note, for ARMA(1,1) we limit the maximum number of past values used in prediction to 10 (for speed)

# Step 1: compute VaR
VaR_T5 <- matrix(NA, nrow = n-1, ncol = length(alpha))
dimnames(VaR_T5) <- list(paste("T", 1:(n-1), sep =""), paste(alpha))
VaR_Joe5 <- VaR_T5
VaR_ast5 <- VaR_T5
VaR_ast11 <- VaR_T5
VaR_clayton11 <- VaR_T5
for(i in 1:(n-1)){
  cat(i,"\n")
  VaR_T5[i,] <- predict(T5fit@tscopula, tail(Ut[1:i], 10), alpha, type = "qf")
  VaR_Joe5[i,] <- predict(joeAR5@tscopula, tail(Ut[1:i], 10), alpha, type = "qf")
  VaR_ast5[i,] <- predict(astAR5@tscopula, tail(Ut[1:i], 10), alpha, type = "qf")
  VaR_clayton11[i,] <- predict(claytonARMA11@tscopula, tail(Ut[1:i], 10), alpha, type = "qf")
  VaR_ast11[i,] <- predict(astARMA11@tscopula, tail(Ut[1:i], 10), alpha, type = "qf")
}

# Step 2: compare visually
ts.plot(qt(Ut,4))
lines(qt(VaR_T5[,1], 4), col = 6)
lines(qt(VaR_Joe5[,1], 4), col = 2)
lines(qt(VaR_ast5[,1], 4), col = 3)
lines(qt(VaR_ast11[,1], 4), col = 4)
lines(qt(VaR_clayton11[,1], 4), col = 5)

# Step 3: compute hit rates
hitrate <- function(v, realized){
  round(100 * mean(realized < v), 2)
}
apply(VaR_T5, 2, FUN = hitrate, realized = Ut[-1])
apply(VaR_Joe5, 2, FUN = hitrate, realized = Ut[-1])
apply(VaR_ast5, 2, FUN = hitrate, realized = Ut[-1])
apply(VaR_ast11, 2, FUN = hitrate, realized = Ut[-1])
apply(VaR_clayton11, 2, FUN = hitrate, realized = Ut[-1])

# Step 4: carry out Christofferson test
# Problem in rugarch code - a large number of hits in a large number of trials can break unconditional 
# test by leading to a NaN
myCCTest <- function(VaRdata, realized, conf.level = 0.99){
  alpha <- as.numeric(dimnames(VaRdata)[[2]])
  output <- rep("", length(alpha))
  names(output) <- paste(alpha)
  for (j in 1: length(alpha))
    output[j] <- myVaRTest(alpha[j], actual = realized, 
                           VaR = VaRdata[,j], 
                           conf.level = conf.level)$uc.Decision
  output
}
myCCTest(VaR_T5, Ut[-1])
myCCTest(VaR_Joe5, Ut[-1])
myCCTest(VaR_ast5, Ut[-1])
myCCTest(VaR_ast11, Ut[-1])
myCCTest(VaR_clayton11, Ut[-1])






