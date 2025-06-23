# Helper function for computing AIC and BIC
logLik.simplefit <- function(output){
  ll <- output$logLik
  attributes(ll) <- list(df = output$npar, nobs = output$nobs)
  class(ll) <- "logLik"
  ll
}

#GARCH/ARCH simulation functions
GARCHvolfunc <- function(x, s, pars, levpar = 0){
  sqrt(pars[1] + (pars[2] + levpar*as.numeric(x < 0))*x^2 + pars[3]*s^2)
}
#
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
fit_tmix <- function(Udata, startpars = c(0.1, 2.5), order = 1, js = FALSE, maxit = 5000){
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


#
# VaRTest function from rugarch (unchanged)
myVaRTest <- function (alpha = 0.05, actual, VaR, conf.level = 0.95) 
{
  N = length(actual)
  VaRn = floor(N * alpha)
  if (N != length(VaR)) 
    stop("\nlength of realized not equal to length of VaR!")
  tmp = LR.cc.test(p = alpha, actual = actual, VaR = VaR, conf.level = conf.level)
  ans = list()
  ans$expected.exceed = floor(alpha * tmp$TN)
  ans$actual.exceed = tmp$N
  ans$uc.H0 = "Correct Exceedances"
  ans$uc.LRstat = tmp$stat.uc
  ans$uc.critical = tmp$crit.val.uc
  ans$uc.LRp = tmp$p.value.uc
  ans$uc.Decision = ifelse(ans$uc.LRp < (1 - conf.level), "Reject H0", 
                           "Fail to Reject H0")
  ans$cc.H0 = "Correct Exceedances & Independent"
  ans$cc.LRstat = tmp$stat.cc
  ans$cc.critical = tmp$crit.val.cc
  ans$cc.LRp = tmp$p.value.cc
  ans$cc.Decision = ifelse(ans$cc.LRp < (1 - conf.level), "Reject H0", 
                           "Fail to Reject H0")
  return(ans)
}

# from rugarch (unchanged)
LR.cc.test <- function (p, actual, VaR, conf.level = 0.95) 
{
  result = .LR.cc(p = p, actual = actual, VaR = VaR)
  crit.val.uc = qchisq(conf.level, df = 1)
  crit.val.cc = qchisq(conf.level, df = 2)
  p.value.cc = 1 - pchisq(result$stat.cc, df = 2)
  p.value.uc = 1 - pchisq(result$stat.uc, df = 1)
  reject = ifelse(p.value.cc < 1 - conf.level, TRUE, FALSE)
  return(list(stat.cc = result$stat.cc, stat.uc = result$stat.uc, 
              p.value.cc = p.value.cc, p.value.uc = p.value.uc, conf.level = conf.level, 
              reject = reject, N = result$N, TN = result$TN, crit.val.uc = crit.val.uc, 
              crit.val.cc = crit.val.cc))
}

# from rugarch (unchanged)
.LR.cc <- function (p, actual, VaR) 
{
  VaR.ind = as.numeric(ifelse(actual < VaR, 1, 0))
  N = sum(VaR.ind)
  TN = length(VaR.ind)
  n <- table(head(VaR.ind, -1), tail(VaR.ind, -1))
  N00 = n[1, 1]
  N11 = n[2, 2]
  N01 = n[1, 2]
  N10 = n[2, 1]
  p01 = N01/(N00 + N01)
  p11 = N11/(N10 + N11)
  pUc = (N01 + N11)/sum(n)
  l1 = (1 - pUc)^(N00 + N10) * pUc^(N01 + N11)
  l2 = (1 - p01)^(N00) * p01^(N01) * (1 - p11)^(N10) * p11^(N11)
  stat.ind = -2 * log(l1/l2)
  stat.uc = .LR.uc(p = p, TN = TN, N = N)
  stat.cc = stat.uc + stat.ind
  return(list(stat.cc = stat.cc, stat.uc = stat.uc, N = N, 
              TN = TN))
}

# This function has been changed as indicated
.LR.uc <- function (p, TN, N) 
{
  #l1 = (1 - p)^(TN - N) * p^N
  # can give rise to NaN
  l1 <- (TN - N) * log(1-p) + N * log(p)
  #l2 = (1 - N/TN)^(TN - N) * (N/TN)^N
  # can give rise to NaA
  l2 <- (TN - N) * log (1 - N/TN) + N * log(N/TN)
  #stat.uc = -2 * log(l1/l2)
  stat.uc <- -2* (l1 - l2)
  return(stat.uc)
}