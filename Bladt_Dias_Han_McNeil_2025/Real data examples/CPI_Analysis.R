library(forecast)
library(tscopula)
#source("simstudyfunctions.R")


# Taken from OECD web page

data <- read.csv(file = "US_CPI_quarterly_raw.csv")
cpivals <- data$Value
dates <- data$TIME
year <- substring(dates,1,4)
qtr <- substring(dates,7,7)
zoodates <- zoo::as.yearqtr(paste(year, qtr, sep="-"))
UScpi_quarterly <- xts::xts(cpivals, zoodates)

# Take data starting in last quarter 1959 (gives data from after 1960 when differenced)
UScpi_quarterly <- UScpi_quarterly[time(UScpi_quarterly) >= "1959 Q4"]
head(UScpi_quarterly)
tail(UScpi_quarterly)
length(UScpi_quarterly)
plot(UScpi_quarterly)

Y <- 400*diff(log(UScpi_quarterly))[-1]
# pdf("alldata.pdf",width = 7, height = 3.5)
plot(Y, main ="")
# dev.off()
Ytrain <- Y[time(Y) <= "2011 Q4"]
length(Ytrain)

# STAGE A

auto.arima(Ytrain, stepwise=FALSE, approximation=FALSE)
Xtrain <- diff(Ytrain)[-1]
(basemod <- auto.arima(Xtrain, stepwise=FALSE, approximation=FALSE))
# Note this gives the same
auto.arima(qnorm(strank(Xtrain)), stepwise=FALSE, approximation=FALSE)

# same model fitted with tscopula for later purposes
# fully Gaussian base model for reference
basemod2 <- tscm(copmod, margin("gauss0"))
(basemod2 <- fit(basemod2, Xtrain, method = "full"))
AIC(basemod, basemod2) 

# STAGE B: Gaussian semiparametric model

U <- strank(Xtrain)
copmod <- sarmacopula(pars = list(ar = 0, ma = c(0, 0), sar = 0))
(modG <- fit(copmod, U))
modG <- fit(modG, U, tsoptions = list(hessian = TRUE))
modG
AIC(modG)
BIC(modG)

# STAGE C: Non-Gaussian semiparametric model

modNG <- auto_dvine(modG,
                         nreplace = 20,
                         tautol = 1e-03,
                         maxlag = Inf, # difference
                         ICtol = 0.5,
                         nstrike = 3,
                         criterion = "AIC",
                         choices = c("gumbel", "clayton", "frank", "joe", "t"),
                         verbose = TRUE)
modNG <- fit(modNG, U, tsoptions = list(hessian = TRUE))
modNG

AIC(modG, modNG)
BIC(modG, modNG)

# Parametric model

margins <- c("norm", "st", "sst", "doubleweibull", "sdoubleweibull",
             "hyp", "gh", "nig", "jsu")
margin_results <- sapply(margins, function(margname, data){
  fit(margin(margname), Xtrain, control = list(maxit = 2000))
}, data = Xdata, USE.NAMES = TRUE)
sapply(margin_results, function(res){res@fit$convergence})
(AIC_margins <- sapply(margin_results, AIC))
bestmargin_index <- which(AIC_margins == min(AIC_margins))
(bestmargin <- margin_results[[bestmargin_index]])
plot(bestmargin)

margmod <- fit(margin("st"), Xtrain)
plot(margmod)
U <- pmarg(bestmargin, Xtrain)
copmod <- sarmacopula(pars = list(ar = 0, ma = c(0, 0), sar = 0))
(copPG <- fit(copPG, U, tsoptions = list(hessian = TRUE)))
AIC(copPG)
BIC(copPG)

copPNG <- auto_dvine(copPG,
                     nreplace = 10,
                     tautol = 1e-03,
                     maxlag = Inf, # difference
                     ICtol = 0.5,
                     nstrike = 3,
                     criterion = "AIC",
                     choices = c("gumbel", "clayton", "frank", "joe", "t"),
                     verbose = TRUE)
copPNG <- fit(copPNG, U, tsoptions = list(hessian = TRUE))
copPNG
AIC(copPNG)
BIC(copPNG)

modPG <- tscm(copPG, margmod)
modPNG <- tscm(copPNG, margmod)

# Summary of copulas

SPtau <- kendall(modNG, 15)
Ptau <- kendall(copPNG, 15)
SPfamnames <- modNG@tscopula@modelspec$family
Pfamnames <- rep("gauss",15)
Pfamnames[1:length(copPNG@tscopula@modelspec$family)] <- copPNG@tscopula@modelspec$family
SPfamnames <- ifelse(SPtau < 0,
                     paste(SPfamnames, modNG@tscopula@modelspec$negrot, sep=""),
                     paste(SPfamnames, modNG@tscopula@modelspec$posrot, sep=""))
Pfamnames <- ifelse(Ptau < 0,
                     paste(Pfamnames, modNG@tscopula@modelspec$negrot, sep=""),
                     paste(Pfamnames, modNG@tscopula@modelspec$posrot, sep=""))

copsummary <- data.frame(SPfamnames,SPtau,Pfamnames,Ptau)
xtable::xtable(copsummary)

# kpacf plot
plot(modNG, plottype = "kendall")
lines((1:30), kendall(copPNG,30), col = "blue")


# forecasting
Xall <- diff(Y)[-1]
Xtest <- Xall[time(Xall) > "2011 Q4"]
Xtrain <- as.numeric(Xtrain)
newdata <- as.numeric(Xtest)

############################

m <- length(newdata)
J <- 20
alpha <- (1:(J-1))/J
SP_NG_scores <- array(NA, dim = c(m, J-1),
                        dimnames = list(forecast = NULL, quantile = paste(alpha)))
SP_G_scores <- SP_NG_scores
P_G_scores <- SP_NG_scores
P_NG_scores <- SP_NG_scores
P_base_scores <- SP_NG_scores
SP_NG_pits <- rep(NA, m)
SP_G_pits <- SP_NG_pits
P_G_pits <- SP_NG_pits
P_NG_pits <- SP_NG_pits
P_base_pits <- SP_NG_pits

xdata <- Xtrain
for (k in 1:length(newdata)){
  xquant_SP_NG <- quantile(xdata, predict(modNG@tscopula, data = strank(xdata), x = alpha, type = "qf"))
  xquant_SP_G <- quantile(xdata, predict(modG@tscopula, data = strank(xdata), x = alpha, type = "qf"))
  xquant_P_G <- predict(modPG, data = xdata, x = alpha, type = "qf")
  xquant_P_NG <- predict(modPNG, data = xdata, x = alpha, type = "qf")
  xquant_P_base <- predict(basemod2, x = alpha, type = "qf")
  SP_NG_scores[k,] <- quantile_score(alpha, xquant_SP_NG, newdata[k])
  SP_G_scores[k,] <- quantile_score(alpha, xquant_SP_G, newdata[k])
  P_NG_scores[k,] <- quantile_score(alpha, xquant_P_NG, newdata[k])
  P_G_scores[k,] <- quantile_score(alpha, xquant_P_G, newdata[k])
  P_base_scores[k,] <- quantile_score(alpha, xquant_P_base, newdata[k])
  #
  SP_NG_pits[k] <- predict(modNG@tscopula, data = strank(xdata), 
                           x = pedf(newdata[k], data = xdata), type = "df")
  SP_G_pits[k] <- predict(modG@tscopula, data = strank(xdata), 
                           x = pedf(newdata[k], data = xdata), type = "df")
  P_G_pits[k] <- predict(basemod2, x = newdata[k], type = "df")
  P_NG_pits[k] <- predict(modPNG, data = xdata, x = newdata[k], type = "df")
  P_G_pits[k] <- predict(modPG, data = xdata, x = newdata[k], type = "df")
  xdata <- c(xdata, newdata[k])
}
# Diebold-Mariano test
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
Diebold_Mariano(SP_G_scores, P_base_scores, dist = "none")
Diebold_Mariano(SP_NG_scores, SP_G_scores, dist = "none")
Diebold_Mariano(P_G_scores, SP_G_scores, dist = "none")
Diebold_Mariano(SP_NG_scores, P_G_scores, dist = "none")
Diebold_Mariano(P_NG_scores, SP_NG_scores, dist = "none")

Diebold_Mariano(SP_G_scores, P_base_scores, dist = "t", df = 3.39)
Diebold_Mariano(SP_NG_scores, SP_G_scores, dist = "t", df = 3.39)
Diebold_Mariano(P_G_scores, SP_G_scores, dist = "t", df = 3.39)
Diebold_Mariano(SP_NG_scores, P_G_scores, dist = "t", df = 3.39)
Diebold_Mariano(P_NG_scores, SP_NG_scores, dist = "t", df = 3.39)


# simple tests
ks.test(SP_NG_pits, "punif")
ks.test(SP_G_pits, "punif")
ks.test(P_G_pits, "punif")
ks.test(P_NG_pits, "punif")
Box.test(P_G_pits, lag = 5, type = "L")
Box.test(P_NG_pits, lag = 5, type = "L")
Box.test(SP_G_pits, lag = 5, type = "L")
Box.test(SP_NG_pits, lag = 5, type = "L")

# Score plot
P_G_ave <- apply(P_G_scores, 2, mean)
P_NG_ave <- apply(P_NG_scores, 2, mean)
SP_G_ave <- apply(SP_G_scores, 2, mean)
SP_NG_ave <- apply(SP_NG_scores, 2, mean)
P_G_ave2 <- P_G_ave/quantile_weight(alpha, "t", df = 3.39)
P_NG_ave2 <- P_NG_ave/quantile_weight(alpha, "t", df = 3.39)
SP_G_ave2 <- SP_G_ave/quantile_weight(alpha, "t", df = 3.39)
SP_NG_ave2 <- SP_NG_ave/quantile_weight(alpha, "t", df = 3.39)

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

# Forecasting plot non-parametric

J <- 20
alpha <- (1:(J-1))/J
NG_forecasts <- array(NA, dim = c(m, J-1),
                   dimnames = list(forecast = NULL, quantile = paste(alpha)))
xdata <- Xtrain
for (k in 1:length(newdata)){
  udata <- strank(xdata)
  uquantNG <- predict(modNG@tscopula, data = udata, x = alpha, type = "qf")
  xquantNG <- quantile(xdata, uquantNG)
  NG_forecasts[k,] <- xquantNG
  xdata <- c(xdata, newdata[k])
}

Ystart <- Y[(time(Y) >= "2011 Q4") & (time(Y) < "2022 Q4")]
Yend <- Y[time(Y) > "2011 Q4"]
Ystartm <- matrix(Ystart, nrow = m, ncol = J-1, byrow= FALSE)
Yforecast <- xts::xts(NG_forecasts + Ystartm, time(Yend))
pp <- plot(Yforecast, col = "firebrick1", lwd = 0.5, main ="")
pp <- xts::addSeries(Yend, on =1, col =1, lwd = 2)
plot(pp)

# Forecasting plot parametric

J <- 20
alpha <- (1:(J-1))/J
P_forecasts <- array(NA, dim = c(m, J-1),
                      dimnames = list(forecast = NULL, quantile = paste(alpha)))
xdata <- Xtrain
for (k in 1:length(newdata)){
  xquantP <- predict(modPNG, data = xdata, x = alpha, type = "qf")
  P_forecasts[k,] <- xquantP
  xdata <- c(xdata, newdata[k])
}

Ystart <- Y[(time(Y) >= "2011 Q4") & (time(Y) < "2022 Q4")]
Yend <- Y[time(Y) > "2011 Q4"]
Ystartm <- matrix(Ystart, nrow = m, ncol = J-1, byrow= FALSE)
Yforecast <- xts::xts(P_forecasts + Ystartm, time(Yend))
pp2 <- plot(Yforecast, col = "blue", lwd = 0.5, main ="")
pp2 <- xts::addSeries(Yend, on =1, col =1, lwd = 2)
plot(pp2)

# Plots for paper

pdf("alldata.pdf",width = 7, height = 3)
plot(Y, main ="")
dev.off()

pdf("Kendall_plot.pdf",width=7, height=4)
par(cex=0.8)
plot(modNG, plottype = "kendall")
lines((1:30), kendall(copPNG,30), col = "blue")
dev.off()

pdf("quantile_score_plot_CPI.pdf",width=7, height=4)
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

pdf("forecasts1.pdf", width =7, height = 3.5)
plot(pp)
dev.off()
pdf("forecasts2.pdf", width =7, height = 4)
plot(pp2)
dev.off()
