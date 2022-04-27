library(tscopula)
library(forecast)

# Preliminary Analysis
data(cpi)
CPI <- 100*diff(log(cpi))[-1]
length(CPI)
plot(CPI)
X <- as.numeric(CPI)
acf(X)
U <- strank(X)

# Automatic arima fit
model <- auto.arima(X, d=0)
model
arcoef <- coef(model)[1:5]
macoef <- coef(model)[6]

# Gauss copula model
modGauss <- armacopula(pars = list(ar = arcoef, ma = macoef))
modGauss <- fit(modGauss, U, tsoptions = list(hessian = TRUE))
modGauss
round(coef(modGauss),3)
round(safe_ses(modGauss@fit$hessian), 3)

# Gumbel copula model
modGumbel <- dvinecopula2(family = "Gumbel",
                          pars = modGauss@tscopula@pars,
                          maxlag = Inf,
                          negtau ="gauss")
modGumbel <- fit(modGumbel, U, tsoptions = list(hessian = TRUE))
modGumbel
round(coef(modGumbel),3)
round(safe_ses(modGumbel@fit$hessian), 3)

# Comparison of residuals
resGauss <- resid(modGauss)
shapiro.test(resGauss)
resGumbel <- resid(modGumbel)
shapiro.test(resGumbel)
plot(modGauss, plottype="residual")
plot(modGumbel, plottype="residual")

# Illustration of kpacf
plot(modGumbel, plottype="kendall")

# Models with margins
modGaussNormal <- tscm(modGauss, margin("norm"))
modGaussNormal <- fit(modGaussNormal, X)
modGaussNormal
modGaussStudent <- tscm(modGauss, margin("sst"))
modGaussStudent <- fit(modGaussStudent, X)
modGaussStudent
modGumbelStudent <- tscm(modGumbel, margin("sst"))
modGumbelStudent <- fit(modGumbelStudent, X)
modGumbelStudent
AIC(modGauss, modGumbel, modGaussNormal, modGaussStudent, modGumbelStudent)
BIC(modGauss, modGumbel, modGaussNormal, modGaussStudent, modGumbelStudent)

# Plot of margin fits
plot(modGaussNormal, plottype="margin")
plot(modGumbelStudent, plottype="margin")
