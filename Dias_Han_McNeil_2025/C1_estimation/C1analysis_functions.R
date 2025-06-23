# Helper function for computing AIC and BIC
logLik.simplefit <- function(output){
  ll <- output$logLik
  attributes(ll) <- list(df = output$npar, nobs = output$nobs)
  class(ll) <- "logLik"
  ll
}

##############################################################################
#GARCH/ARCH simulation functions

GARCHvolfunc <- function(x, s, pars, levpar = 0){
  sqrt(pars[1] + (pars[2] + levpar*as.numeric(x < 0))*x^2 + pars[3]*s^2)
}

simgarch <- function(dist, n, pars, levpar = 0, distpars = NULL, 
                     meanpar = 0, copula = TRUE, outputsigma = FALSE){   
  Z <- switch(dist,
              norm = rnorm(n),
              t = {
                df <- distpars[1]
                rt(n, df)*sqrt((df-2)/df)
              },
              skewt = {
                df <- distpars[1]
                gamma <- distpars[2]
                M1 <- dt(0,df=df)*df/(df-1)*2 
                dist.mn <- M1*(gamma-1/gamma)
                M2 <- df/(df-2)
                dist.var <- (M2-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - M2
                unscaled <- skewt::rskt(n, df = df, gamma= gamma)
                (unscaled-dist.mn)/sqrt(dist.var)
              }
  )
  sigma <- rep(NA, n)
  X <- rep(NA, n)
  sigma[1] <- 1
  X[1] <- Z[1]*sigma[1]
  for (i in 2:n){
    sigma[i] <- GARCHvolfunc(X[i-1], sigma[i-1], pars, levpar) 
    X[i] <- meanpar*X[i-1] + sigma[i] * Z[i]
  }
  if (outputsigma)
    X <- sigma
  if (copula)
    output <- tscopula::strank(X)
  else 
    output <- X
  output
}

######################################################################################
#Estimation helper functions

V2func <- function(u, delta = 0.5, kappa = 1){
  if (delta == 0) 
    return(u)
  else if (delta == 1) 
    return(1 - u)
  else {
    ifelse(u <= delta, 1 - u - (1 - delta) * exp(-kappa * 
                                                   log(delta/u)), u - delta * exp(-(-log((1 - u)/(1 - 
                                                                                                    delta))/kappa)))
  }
}

vtpars <- function(vt, theta, k){
  delta <- switch(vt,
                  symmetric = 0.5,
                  linear = theta[k+1],
                  asymmetric = 0.5,
                  v2p = theta[k+1])
  if ((delta <=0) | (delta >=1))
    delta <- NA
  kappa <- switch(vt,
                  symmetric = 1,
                  linear = 1,
                  asymmetric = theta[k+1],
                  v2p = theta[k+2])
  if (kappa <=0)
    kappa <- NA
  c(delta,kappa)
}

##############################################################################
# Functions for ast copula

dCopula_ast <- function(u, theta){
  dCopula((1+u)/2, tCopula(0, df = theta))
}

hCopula_ast <- function(u, theta, condvar = 1){
  below <- switch(condvar,
                  qt((1+u[,1])/2, df = theta),
                  qt((1+u[,2])/2, df = theta))
  above <- switch(condvar,
                  qt((1+u[,2])/2, df = theta),
                  qt((1+u[,1])/2, df = theta))
  2*pt(sqrt(theta+1)*above/sqrt(theta + below^2), df = theta + 1) -1
}

pCopula_ast <- function(u, theta){
  integrand <- function(u, v, theta){
    hCopula_ast(cbind(u, v), theta)
  }
  mapply(FUN = function(u, v, nu) {
    integrate(integrand,
              lower = 0,
              upper = u,
              v = v,
              theta = nu)$value
  },
  u[,1], u[,2], MoreArgs = list(nu = theta))
}

tauast <- function(nu){
  inner_integrand <- function(v, u, nu){
    pCopula_ast(cbind(u,v), nu) * dCopula_ast(cbind(u,v), nu)
  }
  inner_integral <- function(u, nu){
    sapply(u, function(u, nu){
      integrate(inner_integrand, lower = 0, upper = 1, u = u, nu = nu)$value
    }, nu = nu)}
  4*integrate(inner_integral, lower=0, upper = 1, nu = nu)$value -1
}

dCopula_ast_K <- function(u, theta, alpha1, alpha2){
  powerdata <- cbind(u[,1]^(1-alpha1), u[,2]^(1-alpha2))
  output <- (1-alpha1) * (1-alpha2) * dCopula_ast(powerdata, theta) +
    alpha2 * (1-alpha1) * hCopula_ast(powerdata, theta, condvar = 1)*(u[,2]^(alpha2-1)) +
    alpha1 * (1-alpha2) * hCopula_ast(powerdata, theta, condvar = 2)*(u[,1]^(alpha1-1)) +
    alpha1 * alpha2 * pCopula_ast(powerdata, theta)*(u[,1]^(alpha1-1))*(u[,2]^(alpha2-1))
  output
}

negloglik_iv_ast <- function(theta, data, Khoudraji, v1, v2){
  if (theta[1] <= 0) 
    return(NA)
  indx <- 1
  if (Khoudraji){
    if ((theta[2] < 0) | (theta[2] > 1) | (theta[3] < 0) | (theta[3] > 1))
      return(NA)
    indx <- indx + 2
  }
  vpars1 <- vtpars(v1, theta, indx)
  indx  <- indx + switch(v1,
                         symmetric = 0,
                         linear = 1,
                         asymmetric = 1,
                         v2p = 2)
  vpars2 <- vtpars(v2, theta, indx)
  if (is.na(sum(vpars1 + vpars2)))
    return(NA)
  data[,1] <- V2func(data[,1], delta = vpars1[1], kappa = vpars1[2])
  data[,2] <- V2func(data[,2], delta = vpars2[1], kappa = vpars2[2])
  if (Khoudraji)
    output <- -sum(log(dCopula_ast_K(data, theta[1], theta[2], theta[3])))
  else
    output <- -sum(log(dCopula_ast(data, theta[1])))
  output
}

fit_iv_ast <- function(Udata, Khoudraji = FALSE, v1 = "symmetric", v2 = "symmetric", startpars = NA){
  if (is.na(sum(startpars))){
    startpars <- 2
    if (Khoudraji)
      startpars <- c(startpars, 0.05, 0.05)
    if (v1 != "symmetric")
      startpars <- c(startpars, switch(v1,
                                       asymmetric = 1,
                                       linear = 0.5,
                                       v2p = c(0.5, 1)))
    if (v2 != "symmetric")
      startpars <- c(startpars,switch(v2,
                                      asymmetric = 1,
                                      linear = 0.5,
                                      v2p = c(0.5, 1)))
  }
  data <- cbind(Udata[-length(Udata)],Udata[-1])
  results <- optim(startpars, fn = negloglik_iv_ast, 
                   data = data, Khoudraji = Khoudraji, v1 = v1, v2 = v2, 
                   control = list(warn.1d.NelderMead = FALSE, maxit = 10000))
  par.ests <- results$par
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), 
                 nobs = length(Udata), conv = results$convergence)
  class(output) <- "simplefit"
  output
}


#############################################################################
# Fitting iv-archimedean including Khoudraji

negloglik_iv_archimedean <- function(theta, data, model, 
                                     Khoudraji = FALSE, 
                                     rotate = FALSE, 
                                     v1 = "symmetric",
                                     v2 = "symmetric"){
  #  browser()
  cop <- try(switch(model,
                    Gumbel = copula::gumbelCopula(theta[1]),
                    Clayton = copula::claytonCopula(theta[1]),
                    Joe = copula::joeCopula(theta[1])), silent = TRUE)
  if (inherits(cop, "try-error"))
    return(NA)
  indx <- 1
  if (rotate)
    cop <- copula::rotCopula(cop)
  if (Khoudraji){
    if ((theta[2] < 0) | (theta[2] > 1) | (theta[3] < 0) | (theta[3] > 1))
      return(NA)
    cop <- copula::khoudrajiCopula(copula1 = cop, shapes = c(theta[2], theta[3]))
    indx <- indx + 2
  }
  vpars1 <- vtpars(v1, theta, indx)
  indx  <- indx + switch(v1,
                         symmetric = 0,
                         linear = 1,
                         asymmetric = 1,
                         v2p = 2)
  vpars2 <- vtpars(v2, theta, indx)
  if (is.na(sum(vpars1 + vpars2)))
    return(NA)
  data[,1] <- V2func(data[,1], delta = vpars1[1], kappa = vpars1[2])
  data[,2] <- V2func(data[,2], delta = vpars2[1], kappa = vpars2[2])
  -sum(log(dCopula(data,copula = cop)))
}


fit_iv_archimedean <- function(Udata, model = "Joe", 
                             Khoudraji = FALSE, 
                             rotate = FALSE, 
                             v1 = "symmetric",
                             v2 = "symmetric",
                             startpars = NA){
  if (is.na(sum(startpars))){
    startpars <- 1.2
    if (Khoudraji)
      startpars <- c(startpars, 0.05, 0.05)
    if (v1 != "symmetric")
      startpars <- c(startpars, switch(v1,
                                       asymmetric = 1,
                                       linear = 0.5,
                                       v2p = c(0.5, 1)))
    if (v2 != "symmetric")
      startpars <- c(startpars,switch(v2,
                                      asymmetric = 1,
                                      linear = 0.5,
                                      v2p = c(0.5, 1)))
  }
  data <- cbind(Udata[-length(Udata)],Udata[-1])
  results <- optim(startpars, fn = negloglik_iv_archimedean, 
                   data = data, 
                   model = model,
                   Khoudraji = Khoudraji,
                   rotate = rotate,
                   v1 = v1,
                   v2 = v2,
                   control = list(warn.1d.NelderMead = FALSE, maxit =2000))
  par.ests <- results$par
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), 
       nobs = length(Udata), conv = results$convergence)
  class(output) <- "simplefit"
  output
}

################################################################################
# t mixtures
checkpar_tmix <- function(theta){
  xi1 <- theta[1]
  nu1 <- theta[2]
  xi2 <- theta[3]
  nu2 <- theta[4]
  w <- theta[5]
  ((xi1 <= 0) | (xi1 >= 1) |  (nu1 <= 2) | (nu1 >= 40) |
      (xi2 <= 0) | (xi2 >= 1) |  (nu2 <= 2) | (nu2 >= 40) |
      (w <=0) | (w >= 1))
}
#
dtmix <- function(u, v, theta){
  dcop1 <- rvinecopulib::dbicop(cbind(u, v), family = "t", parameters = c(theta[1], theta[2]))
  dcop2 <- rvinecopulib::dbicop(cbind(u, v), family = "t", parameters = c(-theta[3], theta[4]))
  theta[5] * dcop1 + (1 - theta[5]) * dcop2
}
#
htmix <- function(u, v, theta, cond_var){
  dcop1 <- rvinecopulib::hbicop(cbind(u, v), family = "t", parameters = c(theta[1], theta[2]), cond_var)
  dcop2 <- rvinecopulib::hbicop(cbind(u, v), family = "t", parameters = c(-theta[3], theta[4]), cond_var)
  theta[5] * dcop1 + (1 - theta[5]) * dcop2
}
#
negloglik_tmix <- function(theta, data, order, js){
  loglik <- 0
  for (i in 1:order){
    if (js){
      thetai <- (i-1)*2 + 1:2
      thetavals <- c(theta[thetai], theta[thetai], 0.5)
    }
    else{
      thetai <- (i-1)*5 + 1:5
      thetavals <- theta[thetai]
    }
    if (checkpar_tmix(thetavals))
      return(NA)
    loglik <- loglik + sum(log(dtmix(data[,1], data[,2], thetavals)))
    if (i < order){
      u1 <- htmix(data[,1], data[,2], thetavals, 2)
      u2 <- htmix(data[,1], data[,2], thetavals, 1)
      data <- cbind(u1[-length(u1)],u2[-1])
    }
  }
  -loglik
}
#
fit_tmix <- function(Udata, startpars, order = 1, js = FALSE, maxit = 5000){
  if ((is.na(sum(startpars))) & js)
    startpars <- c(0.3, 2.4)
  n <- length(Udata)
  data <- cbind(Udata[-n], Udata[-1])
  results <- optim(startpars, fn = negloglik_tmix, 
                   data = data, order = order, js = js, control = list(maxit = maxit))
  par.ests <- results$par
  conv <- results$convergence
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), nobs = n, conv = conv)
  class(output) <- "simplefit"
  output
}
#
################################################################################
# Archimedean mixtures
# Fitting jointly symmetric Archimedean mixture

negloglik_mix_archimedean <- function(theta, data, model, Khoudraji, rotate){
  cop <- try(switch(model,
                    Gumbel = copula::gumbelCopula(theta[1]),
                    Clayton = copula::claytonCopula(theta[1]),
                    Joe = copula::joeCopula(theta[1])), silent = TRUE)
  if (inherits(cop, "try-error"))
    return(NA)
  if (Khoudraji){
    if ((theta[2] < 0) | (theta[2] > 1) | (theta[3] < 0) | (theta[3] > 1))
      return(NA)
    cop <- copula::khoudrajiCopula(copula1 = cop, shapes = c(theta[2], theta[3]))
  }
  -sum(log(
    (
      0.25*copula::dCopula(data, cop)
      +0.25*copula::dCopula(data, copula::rotCopula(cop, flip = c(TRUE, FALSE)))
      +0.25*copula::dCopula(data, copula::rotCopula(cop, flip = c(FALSE, TRUE)))
      +0.25*copula::dCopula(data, copula::rotCopula(cop, flip = c(TRUE, TRUE)))
    )
  ))
}
#
fit_mix_archimedean <- function(Udata, model, Khoudraji = FALSE, startpars = NA){
  if (is.na(sum(startpars))){
    startpars <- c(1.2)
    if (Khoudraji) 
      startpars <- c(startpars, 0.5, 0.5)
  }
  data <- cbind(Udata[-length(Udata)],Udata[-1])
  results <- optim(startpars, fn = negloglik_mix_archimedean, 
                   data = data, 
                   model = model,
                   Khoudraji = Khoudraji,
                   control = list(warn.1d.NelderMead = FALSE))
  par.ests <- results$par
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), 
                 nobs = length(Udata), conv = results$convergence)
  class(output) <- "simplefit"
  output
}
###############################################################################
# convex Gumbel functions 
checkpar_CGmix <- function(theta){
  w1 <- theta[1]
  par1 <- theta[3]
  w2 <- theta[2]
  par2 <- theta[4]
  w <- theta[5]
  ((w1 <= 0) | (w1 >= 1) |  (par1 <= 1) |
      (w2 <= 0) | (w2 >= 1) |  (par2 <= 1) |
      (w <=0) | (w >= 1))
}
#
dCGmix <- function(u, v, theta){
  data <- cbind(u,v)
  theta[5] * theta[1] * rvinecopulib::dbicop(data, family = "gumbel", parameters = theta[3]) +
    theta[5]* (1-theta[1]) * rvinecopulib::dbicop(data, family = "gumbel", rotation = 180, parameters = theta[3]) +
    (1 - theta[5]) * theta[2] * rvinecopulib::dbicop(data, family = "gumbel", rotation = 90, parameters = theta[4]) +
    (1 - theta[5]) * (1-theta[2]) * rvinecopulib::dbicop(data, family = "gumbel", rotation = 270, parameters = theta[4])
}
#
hCGmix <- function(u, v, theta, cond_var){
  data <- cbind(u,v)
  theta[5] * theta[1] * rvinecopulib::hbicop(data, family = "gumbel", parameters = theta[3], cond_var) +
    theta[5]* (1-theta[1]) * rvinecopulib::hbicop(data, family = "gumbel", rotation = 180, parameters = theta[3], cond_var) +
    (1 - theta[5]) * theta[2] * rvinecopulib::hbicop(data, family = "gumbel", rotation = 90, parameters = theta[4], cond_var) +
    (1 - theta[5]) * (1-theta[2]) * rvinecopulib::hbicop(data, family = "gumbel", rotation = 270, parameters = theta[4], cond_var)
}
negloglik_CGmix <- function(theta, data, order){
  loglik <- 0
  for (i in 1:order){
    thetai <- (i-1)*5 + 1:5
    thetavals <- theta[thetai]
    if (checkpar_CGmix(thetavals))
      return(NA)
    loglik <- loglik + sum(log(dCGmix(data[,1], data[,2], thetavals)))
    if (i < order){
      u1 <- hCGmix(data[,1], data[,2], thetavals, 2)
      u2 <- hCGmix(data[,1], data[,2], thetavals, 1)
      data <- cbind(u1[-length(u1)],u2[-1])
    }
  }
  -loglik
}
#
fit_CGmix <- function(Udata, startpars, order = 1, maxit = 5000){
  data <- cbind(Udata[-length(Udata)],Udata[-1])
  results <- optim(startpars, fn = negloglik_CGmix, 
                   data = data, order = order, control = list(maxit = maxit))
  par.ests <- results$par
  n <- dim(data)[1]+1
  conv <- results$convergence
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), nobs = n, conv = conv)
  class(output) <- "simplefit"
  output
}
############################################################################
# Bespoke functions

fit_mix_t <- function(data, startpars = NA){
  fit_tmix(data, startpars = startpars, order = 1, js = TRUE)
}

fit_mix_Clayton <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Clayton", startpars = startpars)
}

fit_mix_Joe <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Joe", startpars = startpars)
}

fit_mix_Gumbel <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Gumbel", startpars = startpars)
}

fit_iv_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", rotate = TRUE, startpars = startpars)
}

fit_iv_Gumbel <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Gumbel", startpars = startpars)
}

fit_iv_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", startpars = startpars)
}

fit_mix_K_Clayton <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Clayton", Khoudraji = TRUE, startpars = startpars)
}

fit_mix_K_Joe <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Joe", Khoudraji = TRUE, startpars = startpars)
}

fit_mix_K_Gumbel <- function(data, startpars = NA){
  fit_mix_archimedean(data, model = "Gumbel", Khoudraji = TRUE, startpars = startpars)
}

fit_iv_K_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", rotate = TRUE, Khoudraji = TRUE, startpars = startpars)
}

fit_iv_K_Gumbel <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Gumbel", Khoudraji = TRUE, startpars = startpars)
}

fit_iv_K_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", Khoudraji = TRUE, startpars = startpars)
}

fit_iv_K_ast <- function(data, startpars = NA){
  fit_iv_ast(data, Khoudraji = TRUE, startpars = startpars)
}

fit_iva_ast <- function(data, startpars = NA){
  fit_iv_ast(data, v1 = "asymmetric", startpars = startpars)
}

fit_iva_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", rotate = TRUE, v1 = "asymmetric", startpars = startpars)
}

fit_iva_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", v1 = "asymmetric", startpars = startpars)
}

fit_iva_K_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", 
                     rotate = TRUE, Khoudraji = TRUE, v1 = "asymmetric", startpars = startpars)
}

fit_iva_K_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", Khoudraji = TRUE, v1 = "asymmetric", startpars = startpars)
}

fit_iva_K_ast <- function(data, startpars = NA){
  fit_iv_ast(data, Khoudraji = TRUE, v1 = "asymmetric", startpars = startpars)
}

fit_ivl_ivl_ast <- function(data, startpars = NA){
  fit_iv_ast(data, v1 = "linear", v2 = "linear", startpars = startpars)
}

fit_ivl_ivl_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", rotate = TRUE, v1 = "linear", v2="linear", startpars = startpars)
}

fit_ivl_ivl_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", v1 = "linear", v2="linear", startpars = startpars)
}

fit_ivl_ivl_K_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", 
                     rotate = TRUE, Khoudraji = TRUE, v1 = "linear", v2="linear", startpars = startpars)
}

fit_ivl_ivl_K_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", Khoudraji = TRUE, v1 = "linear", v2="linear", startpars = startpars)
}

fit_ivl_ivl_K_ast <- function(data, startpars = NA){
  fit_iv_ast(data, Khoudraji = TRUE, v1 = "linear", v2="linear", startpars = startpars)
}

fit_iv2_iv2_ast <- function(data, startpars = NA){
  fit_iv_ast(data, v1 = "v2p", v2 = "v2p", startpars = startpars)
}

fit_iv2_iv2_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", rotate = TRUE, v1 = "v2p", v2="v2p", startpars = startpars)
}

fit_iv2_iv2_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", v1 = "v2p", v2="v2p", startpars = startpars)
}

fit_iv2_iv2_K_Clayton_180 <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Clayton", 
                     rotate = TRUE, Khoudraji = TRUE, v1 = "v2p", v2="v2p", startpars = startpars)
}

fit_iv2_iv2_K_Joe <- function(data, startpars = NA){
  fit_iv_archimedean(data, model = "Joe", Khoudraji = TRUE, v1 = "v2p", v2="v2p", startpars = startpars)
}

fit_iv2_iv2_K_ast <- function(data, startpars = NA){
  fit_iv_ast(data, Khoudraji = TRUE, v1 = "v2p", v2="v2p", startpars = startpars)
}





