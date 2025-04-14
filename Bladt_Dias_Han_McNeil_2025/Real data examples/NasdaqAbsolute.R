library(forecast)
library(tscopula)
#source("simstudyfunctions.R")

# compute quantile scores
quantile_score <- function(alpha, quants, forecast){
  weights <- ifelse(forecast <= quants, 1-alpha, alpha)
  weights*abs(forecast - quants)
}


quantile_weight <- function(theta, dist = "norm", inv = FALSE, ...){
  n <- length(theta)
  output <- rep(NA, n)
  integrand <- function(x, alpha, dist, ...){
    qfunc <- eval(parse(text = paste("q", dist, sep ="")))
    (as.numeric(qfunc(x, ...) <= qfunc(alpha, ...))-alpha)*(qfunc(alpha, ...) - qfunc(x, ...))
  }
  for (i in 1:n)
    output[i] <- integrate(integrand, lower = 0, upper = 1, 
                           alpha = theta[i], dist = dist, ...)$value
  if (inv)
    output <- 1/output
  output
}

Diebold_Mariano <- function(scores1, scores2, dist = "unif", ...){
  if (!(identical(dim(scores1), dim(scores2))))
    stop ("Score matrices of unequal size")
  m <- dim(scores1)[1]
  alpha <- as.numeric(dimnames(scores1)[[2]])
  wts <- switch(dist,
                "none" = rep(1,length(alpha)),
                "unif" = 1/(alpha*(1-alpha)/2),
                1/quantile_weight(alpha, dist, ...))
  wtm <- matrix(wts, nrow = m, ncol = length(alpha), byrow=TRUE) 
  AWQS1i <- apply(scores1*wtm,1, mean)
  AWQS2i <- apply(scores2*wtm,1, mean)
  sigma <- sqrt(mean((AWQS1i - AWQS2i)^2))
  AWQS1 <- mean(AWQS1i)
  AWQS2 <- mean(AWQS2i)
  Tstat <- sqrt(m)*(AWQS1 - AWQS2)/sigma
  pvalue <- 2*(1-pnorm(abs(Tstat)))
  bestmethod <- paste("Method",ifelse(Tstat <= 0, 1, 2))
  list(AWQS1 = AWQS1, AWQS2 = AWQS2, sigma = sigma, Tstat=Tstat, Better = bestmethod, Null = "Forecast methods equally accurate",
       pvalue=pvalue)
}

# quantmod::getSymbols("^IXIC", src = "yahoo", auto.assign = TRUE)
# save("nasdaq", file = "nasdaq.RData")
load("nasdaq.RData")
plot(nasdaq)
X <- abs((diff(log(nasdaq))[-1]))*100
X <- X[X!=0]
plot(X)
Xtrain <- X["/2019-12-31"] # pre-covid period
Xtest <- X["2020-01-01/2020-12-31"] # covid and the stock market crash
length(Xtrain)
length(Xtest)

pdf("nasdaqdata.pdf",width = 7, height = 3)
plot(X["/2020-12-31"], main ="")
dev.off()

# Stage A

U <- strank(Xtrain)
(basemod <- auto.arima(qnorm(U), stationary = TRUE, stepwise=FALSE, approximation=FALSE))

# STAGE B: Gaussian semiparametric model

copmod <- armacopula(pars = list(ar = c(0), ma = c(0)))
(cop_SP_G <- fit(copmod, U, tsoptions = list(hessian = TRUE), control = list(maxit = 2000)))

# COMPARE with parametric

margmod <- fit(new("margin", name = "gamma", pars = c(shape=1, scale =1)), Xtrain)
plot(margmod)
margmod2 <- fit(new("margin", name = "weibull", pars = c(shape=1, scale =1)), Xtrain)

dig <- function(x, chi, psi, logvalue){
  dgig(x, lambda = -0.5, abs(chi), abs(psi), logvalue)
}
dinvgamma <- function(x, shape, rate, logvalue){
  dgig(x, lambda = -abs(shape), chi = 2*abs(rate), psi=0, logvalue)
}
margmod3 <- fit(new("margin", name = "ig", pars = c(chi= 0.5, psi= 0.5)), Xtrain)
margmod4 <- fit(new("margin", name = "invgamma", pars = c(shape=1, rate =1)), Xtrain)
AIC(margmod, margmod2, margmod3,margmod4)
BIC(margmod, margmod2, margmod3,margmod4)

plot(margmod)
UP <- pmarg(margmod, Xtrain)
(cop_P_G <- fit(copmod, UP, tsoptions = list(hessian = TRUE)))

pdf("QQplot.pdf",width=4, height=4)
par(cex=0.8, pty = "s")
plot(margmod)
dev.off()

# STAGE C: Non-Gaussian semiparametric model

cop_SP_NG <- auto_dvine(cop_SP_G,
                    nreplace = 10,
                    tautol = 1e-03,
                    maxlag = Inf, # difference
                    ICtol = 0.2,
                    nstrike = 3,
                    criterion = "AIC",
                    choices = c("gumbel", "clayton", "frank", "joe", "t"),
                    verbose = TRUE)
cop_SP_NG <- fit(cop_SP_NG, U, tsoptions = list(hessian = TRUE))
cop_SP_NG

# Non-Gaussian parametric model

cop_P_NG <- auto_dvine(cop_P_G,
                     nreplace = 10,
                     tautol = 1e-03,
                     maxlag = Inf, # difference
                     ICtol = 0.2,
                     nstrike = 3,
                     criterion = "AIC",
                     choices = c("gumbel", "clayton", "frank", "joe", "t"),
                     verbose = TRUE)
cop_P_NG <- fit(cop_P_NG, UP, tsoptions = list(hessian = TRUE))
cop_P_NG

AIC(cop_P_G, cop_SP_G, cop_SP_NG, cop_P_NG)
BIC(cop_P_G, cop_SP_G, cop_SP_NG, cop_P_NG)

# Residual plots
plot(cop_SP_NG)
plot(cop_P_NG)

# Residual test
res_SP_NG<- resid(cop_SP_NG)
shapiro.test(res_SP_NG)
res_P_NG<- resid(cop_P_NG)
shapiro.test(res_P_NG)

# Make full parametric models for forecasting
mod_P_G <- tscm(cop_P_G, margmod)
mod_P_NG <- tscm(cop_P_NG, margmod)

# Forecasting
newdata <- as.numeric(Xtest)

m <- length(newdata)
J <- 20
alpha <- (1:(J-1))/J
P_G_scores <- array(NA, dim = c(m, J-1),
                   dimnames = list(forecast = NULL, quantile = paste(alpha)))
P_NG_scores <- P_G_scores
SP_G_scores <- P_G_scores
SP_NG_scores <- P_G_scores


xdata <- as.numeric(Xtrain)
for (k in 1:length(newdata)){
  xquant_SP_NG <- quantile(xdata, predict(cop_SP_NG@tscopula, data = strank(xdata), x = alpha, type = "qf"))
  xquant_SP_G <- quantile(xdata, predict(cop_SP_G@tscopula, data = strank(xdata), x = alpha, type = "qf"))
  xquant_P_G <- predict(mod_P_G, data = xdata, x = alpha, type = "qf")
  xquant_P_NG <- predict(mod_P_NG, data = xdata, x = alpha, type = "qf")
  SP_NG_scores[k,] <- quantile_score(alpha, xquant_SP_NG, newdata[k])
  SP_G_scores[k,] <- quantile_score(alpha, xquant_SP_G, newdata[k])
  P_NG_scores[k,] <- quantile_score(alpha, xquant_P_NG, newdata[k])
  P_G_scores[k,] <- quantile_score(alpha, xquant_P_G, newdata[k])
  xdata <- c(xdata, newdata[k])
}

# Diebold-Mariano tests
Diebold_Mariano(P_G_scores, SP_G_scores, dist = "none")
Diebold_Mariano(SP_NG_scores, P_G_scores, dist = "none")
Diebold_Mariano(P_NG_scores, SP_NG_scores, dist = "none")

Diebold_Mariano(P_G_scores, SP_G_scores, dist = "gamma", shape = 1, scale = 0.7)
Diebold_Mariano(SP_NG_scores, P_G_scores, dist = "gamma", shape = 1, scale = 0.7)
Diebold_Mariano(P_NG_scores, SP_NG_scores, dist = "gamma", shape = 1, scale = 0.7)

# Order of models - best to worst
# P_NG, SP_NG, P_G, P_NG
# Significant difference between middle 2

# Score plot
P_G_ave <- apply(P_G_scores, 2, mean)
P_NG_ave <- apply(P_NG_scores, 2, mean)
SP_G_ave <- apply(SP_G_scores, 2, mean)
SP_NG_ave <- apply(SP_NG_scores, 2, mean)
P_G_ave2 <- P_G_ave/quantile_weight(alpha, "gamma", shape = 1, scale = 0.7)
P_NG_ave2 <- P_NG_ave/quantile_weight(alpha, "gamma", shape = 1, scale = 0.7)
SP_G_ave2 <- SP_G_ave/quantile_weight(alpha, "gamma", shape = 1, scale = 0.7)
SP_NG_ave2 <- SP_NG_ave/quantile_weight(alpha, "gamma", shape = 1, scale = 0.7)

pdf("quantile_score_plot_nasdaq.pdf",width=7, height=4)
par(mfrow=c(1,2), cex =0.8)
plot(alpha, P_G_ave, type ="l",ylim = range(P_G_ave, P_NG_ave, SP_G_ave, SP_NG_ave), ylab = "quantile score")
lines(alpha, SP_G_ave, col = "black", lty=2)
lines(alpha, SP_NG_ave, col = "red")
lines(alpha, P_NG_ave, col = "blue")
plot(alpha, P_G_ave2, type ="l",ylim = range(P_G_ave2, P_NG_ave2, SP_G_ave2, SP_NG_ave2), ylab = "weighted quantile score")
lines(alpha, SP_G_ave2, col = "black", lty=2)
lines(alpha, SP_NG_ave2, col = "red")
lines(alpha, P_NG_ave2, col = "blue")
par(mfrow=c(1,1))
dev.off()

(SPnms <- cop_SP_NG@tscopula@modelspec$family)
(Pnms <- cop_P_NG@tscopula@modelspec$family)
SPtau <- kendall(cop_SP_NG, length(SPnms))
Ptau <- kendall(cop_P_NG, length(Pnms))
SPnms <- ifelse(SPtau < 0,
                     paste(SPnms, cop_SP_NG@tscopula@modelspec$negrot, sep=""),
                     paste(SPnms, cop_SP_NG@tscopula@modelspec$posrot, sep=""))
Pnms <- ifelse(Ptau < 0,
                paste(Pnms, cop_P_NG@tscopula@modelspec$negrot, sep=""),
                paste(Pnms, cop_P_NG@tscopula@modelspec$posrot, sep=""))

nl <- max(length(SPnms), length(Pnms))
Pfamnames <- rep("gauss", nl)
SPfamnames <- rep("gauss", nl)
Ptau <- kendall(cop_P_NG, nl)
SPtau <- kendall(cop_SP_NG, nl)
Pfamnames[1:length(Pnms)] <- Pnms
SPfamnames[1:length(SPnms)] <- SPnms
(copsummary <- data.frame(SPfamnames,SPtau,Pfamnames,Ptau))
xtable::xtable(copsummary)

plot(cop_P_NG, plottype="kendall")

