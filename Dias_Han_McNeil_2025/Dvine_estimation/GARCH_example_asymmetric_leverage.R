library(tscopula)
library(rvinecopulib)
source("Dvine_functions.R")

#tsoptions <- list(hessian = TRUE, method = "Nelder-Mead", avoidzero= FALSE)
nsim <- 10000

# asymmetric innovations and leverage
nu <- 4
lambda <- 0.8
lev <- 0.4
GARCHpars = c(0.1, 0.1, 0.6)
set.seed(123)
Ut <- simgarch("skewt", n = nsim, pars = GARCHpars, levpar = lev, distpars = c(nu, lambda))

ts.plot(Ut)
# Should be approx uniform
hist(Ut)

# Copula T1-T5
# Fitted with tscopula package
T1fit <- fit(dvinecopula(family = "t", pars = list(c(0, 10))), Ut)
T2fit <- fit(dvinecopula(family = "t", pars = list(c(0, 7), c(0, 10))), Ut)
# Exact results for T5 are sensitive to starting values; slowish fit
T5fit <- fit(dvinecopula(family = "t", 
                       pars = list(c(0, 6),c(0, 7), c(0, 8), c(0, 9), c(0, 10))), Ut, 
           control = list(maxit = 5000))


# Copula B1-B5
B1fit <- fit_CGmix(Ut, startpars = c(0.5, 0.5, 2, 2, 0.5))
B2fit <- fit_CGmix(Ut, startpars = rep(c(0.3, 0.75, 1.2, 1.2, 0.5), 2), order = 2)
# very slow fit
thetaB0 <- c(0.4018446,0.8751692,1.449143,1.548279,0.5083089,
             0.5400843,0.8332324,1.282219,1.251298,0.4480213,
             0.2994759,0.8315824,1.323381,1.091765,0.2487844,
             0.4272739,0.907188,1.221283,1.082203,0.2665131,
             0.3087835,0.7548283,1.120064,1.072388,0.3674199)
B5fit <- fit_CGmix(Ut, startpars = thetaB0, order = 5, maxit = 20000)


# Using inverse-v-transformed copulas

# AR(1)-type models
startpars <- list(ar = 0.3)
joeAR1 <- fit(dvinecopulavt(family = "joe", pars = startpars), Ut)
claytonAR1 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars), Ut)
astAR1 <- fit(dvinecopulavt(family = "ast", pars = startpars), Ut)

# AR(2)-type models
startpars <- list(ar = c(0.25, 0.1))
joeAR2 <- fit(dvinecopulavt(family = "joe", pars = startpars), Ut)
claytonAR2 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars), Ut)
astAR2 <- fit(dvinecopulavt(family = "ast", pars = startpars), Ut)

# AR(5)-type models
startpars <- list(ar = c(0.25, 0.1, 0.05, 0.05, 0.04))
# Joe and Clayton need increased iterations
joeAR5 <- fit(dvinecopulavt(family = "joe", pars = startpars), Ut,
              control= list(maxit = 10000))
claytonAR5 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars), Ut,
                  control= list(maxit = 10000))
astAR5 <- fit(dvinecopulavt(family = "ast", pars = startpars), Ut)

# ARMA(1,1)-type models
startpars <- list(ar = 0.95, ma = -0.85)
joeARMA11 <- fit(dvinecopulavt(family = "joe", pars = startpars, 
                               tautol = 1e-06, maxlag = 20), Ut)
claytonARMA11 <- fit(dvinecopulavt(family = "clayton", rotation = 180, pars = startpars,
                                   tautol = 1e-05, maxlag = 20), Ut)
astARMA11 <- fit(dvinecopulavt(family = "ast", pars = startpars,
                               tautol = 5e-06, maxlag = 20), Ut,
                 tsoptions= list(hessian = TRUE))

# Results
allmods <- list("T1" = T1fit, "B1" = B1fit, 
                "Joe-AR(1)" = joeAR1, "Clayton180-AR(1)" = claytonAR1, "ast-AR(1)" = astAR1,
                "T2" = T2fit, "B2" = B2fit,
                "Joe-AR(2)" = joeAR2, "Clayton180-AR(2)" = claytonAR2, "ast-AR(2)" = astAR2,
                "T5" = T5fit, "B5" = B5fit,
                "Joe-AR(5)" = joeAR5, "Clayton180-AR(5)" = claytonAR5, "ast-AR(5)" = astAR5,
                "Joe-ARMA(1,1" = joeARMA11, "Clayton180-ARMA(1,1)" = claytonARMA11, 
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