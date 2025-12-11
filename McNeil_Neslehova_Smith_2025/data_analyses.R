library(basiscor)
library(copula)
library(corrplot)


# Estimation functions

vt <- function(u, delta, kappa){
  output <- u
  lower <- u[u<= delta]
  upper <- u[u> delta]
  output[u <= delta] <- (1-lower) - (1-delta)*(1- (1- lower/delta)^kappa)
  output[u > delta] <- upper - delta*(1- (1- (1-upper)/(1-delta))^(1/kappa))
  output
}

copdens <- function(data, model, theta, case){
  cop <- try(switch(model,
                    Clayton = copula::rotCopula(copula::claytonCopula(theta[1])),
                    Gumbel = copula::gumbelCopula(theta[1]),
                    Gauss = copula::normalCopula(theta[1]),
                    Frank = copula::frankCopula(theta[1])), silent = TRUE)
  if (inherits(cop, "try-error"))
    return(NA)
  if ((theta[2] <=0) | theta[2] >= 1)
    return(NA)
  if (theta[3] <= 0)
    return(NA)
  if (case > 1){
    if ((theta[4] <=0) | theta[4] >= 1)
      return(NA)
    if (theta[5] <= 0)
      return(NA)
  }
  data[,1] <- vt(data[,1], delta = theta[2], kappa = theta[3])
  if (case == 2) 
    data[,2] <- vt(data[,2], delta = theta[4], kappa = theta[5])
  if (case == 3) 
    data[,2] <- 1 - vt(data[,2], delta = theta[4], kappa = theta[5])
  dCopula(data,copula = cop)
}

negloglik_model <- function(theta, data, model, case){
  -sum(log(copdens(data, model, theta, case)))
}

fit_copula <- function(data, model, case)
{
  if (model == "Gauss")
    startpars <- c(0.5, 0.5, 1)
  else
    startpars <- c(1.2, 0.5, 1)
  if (case > 1)
    startpars <- c(startpars, 0.5, 1)
  results <- optim(startpars, fn = negloglik_model, 
                   data = data, 
                   model = model,
                   case = case,
                   control = list(warn.1d.NelderMead = FALSE, maxit =2000))
  par.ests <- results$par
  output <- list(logLik = -results$value, par.ests = par.ests, npar = length(par.ests), 
                 nobs = dim(data)[1], conv = results$convergence)
  class(output) <- "simplefit"
  output
}

######################################################
# Data
n <- 200
set.seed(15)

X <- rnorm(n)
Y <- X^2 + rnorm(n)
data1 <- cbind(X, Y)
basiscor(data)

X <- rnorm(n)
Y2 <- X^2 + abs(rnorm(n,0,1.5))
sign <- 2*rbinom(n, 1, 0.5) -1
Y <- sign*Y2
data2 <- cbind(X, Y)
basiscor(data2)

theta <- runif(n)*2*pi
eps <- rnorm(length(theta), 0, 0.1)
X <- cos(theta)*(1 + eps)
Y <- sin(theta)*(1 + eps)
data3 <- cbind(X, Y)
basiscor(data3)

pdf("nonmonpic.pdf", width = 7, height = 3)
opar <- par(pty = "s", mar = c(5.1, 4.1, 4.1, 2.1) - 1, 
            mgp=c(0.5,0,0))
lay <- matrix(1:3, ncol = 3, byrow = TRUE) # layout matrix
layout(lay)
plot(data1, axes=FALSE)
box()
plot(data2, axes=FALSE)
box()
plot(data3, axes=FALSE)
box()
dev.off()

###############################################
pdf("chessboards.pdf", width = 7, height = 3)
opar <- par(pty = "s", 
            mar = c(5.1, 4.1, 4.1, 2.1) - 1, 
            mgp=c(1.4,0.5,0))
lay <- matrix(1:3, ncol = 3, byrow = TRUE) # layout matrix
layout(lay)
maxorder <- 8
type <- "legendre"
mat1 <- basiscormatrix(data1, maxorder = maxorder, type = type, method = "T1")
mat2 <- basiscormatrix(data2, maxorder = maxorder, type = type, method = "T1")
mat3 <- basiscormatrix(data3, maxorder = maxorder, type = type, method = "T1")
corrplot(mat1, method = "color")
corrplot(mat2, method = "color")
corrplot(mat3, method = "color")
dev.off()

#########################################
# Estimation
expansion1 <- basisexpand(data1)
expansion2 <- basisexpand(data2)
expansion3 <- basisexpand(data3)
Udata1 <- apply(data1, 2, copula::pobs)
Udata2 <- apply(data2, 2, copula::pobs)
Udata3 <- apply(data3, 2, copula::pobs)
Udata1B <- udpbexdata(Udata1, expansion1, boundaryadjust = TRUE)
# hist(Udata1B[,1])
# hist(Udata1B[,2])
Udata2B <- udpbexdata(Udata2, expansion2, boundaryadjust = TRUE)
Udata3B <- udpbexdata(Udata3, expansion3, boundaryadjust = TRUE)
(copmod1 <- fitCopula(rotCopula(claytonCopula()), Udata1B))
fitCopula(gumbelCopula(), Udata1B)
(copmod2 <- fitCopula(rotCopula(claytonCopula()), Udata2B))
fitCopula(gumbelCopula(), Udata2B)
(copmod3 <- fitCopula(frankCopula(), Udata3B))
fitCopula(normalCopula(), Udata3B)
(parresults1 <- fit_copula(Udata1, "Clayton", 1))
fit_copula(Udata1, "Gumbel", 1)
(parresults2 <- fit_copula(Udata2, "Gumbel", 2))
fit_copula(Udata2, "Clayton", 2)
(parresults3 <- fit_copula(Udata3, "Frank", 3))
fit_copula(Udata3, "Gauss", 3)
#############################################

# small graphical illustrations
par(mfrow=c(1,2))
basiscor:::pcsm_analysis(expansion1$alphag, expansion1$type, plotfunc="TRUE")
abline(h=-0.55)
rts <- basiscor:::root_analysis(-0.55, expansion1$alphag, expansion1$type)
abline(v = rts$roots, lty = 3)
sum(rts$roots * rts$signs)
plot(expansion1,udp = TRUE)
abline(v = rts$roots, lty = 3)
abline(h = sum(rts$roots * rts$signs), lty=3)
par(mfrow=c(1,1))
# small graphical illustration B
par(mfrow=c(1,2))
basiscor:::pcsm_analysis(expansion3$alphah, expansion1$type, plotfunc="TRUE")
abline(h=-0.55)
rts <- basiscor:::root_analysis(-0.55, expansion3$alphah, expansion3$type)
abline(v = rts$roots, lty = 3)
sum(rts$roots * rts$signs) # add one necessary
plot(expansion3,udp = TRUE)
abline(v = rts$roots, lty = 3)
abline(h = sum(rts$roots * rts$signs) + 1, lty=3)
par(mfrow=c(1,1))

###############################################
pdf("basisexpansions.pdf", width = 7, height = 3)
opar <- par(pty = "s", 
            mar = c(5.1, 4.1, 4.1, 2.1) - 1, 
            mgp=c(1.4,0.5,0))
lay <- matrix(1:3, ncol = 3, byrow = TRUE) # layout matrix
layout(lay)
plot(expansion1)
plot(expansion2)
plot(expansion3)
dev.off()



pdf("estimationdata.pdf", width = 7, height = 7.5)
opar <- par(pty = "s", 
            mar = c(5.1, 4.1, 4.1, 2.1) - 2, 
            mgp=c(1.4,0.5,0))
lay <- matrix(1:9, ncol = 3, byrow = TRUE) # layout matrix
layout(lay)
plot(Udata1, xlab = "U", ylab = "V")
plot(Udata2, xlab = "U", ylab = "V")
plot(Udata3, xlab = "U", ylab = "V")
plot(expansion1, udp = TRUE)
u <- seq(from=0, to = 1, length=200)
lines(u, vt(u, parresults1$par.ests[2], parresults1$par.ests[3]), col = "black", lty=3)
lines(u, u, col = "red", lty=3)
plot(expansion2, udp = TRUE)
lines(u, vt(u, parresults2$par.ests[2], parresults2$par.ests[3]), col = "black", lty=3)
lines(u, vt(u, parresults2$par.ests[4], parresults2$par.ests[5]), col = "red", lty=3)
plot(expansion3, udp = TRUE)
lines(u, vt(u, parresults3$par.ests[2], parresults3$par.ests[3]), col = "black", lty=3)
lines(u, 1-vt(u, parresults3$par.ests[4], parresults3$par.ests[5]), col = "red", lty=3)
plot(Udata1B, xlab = expression(T[g](U)), ylab = expression(T[h](V)))
plot(Udata2B, xlab = expression(T[g](U)), ylab = expression(T[h](V)))
plot(Udata3B, xlab = expression(T[g](U)), ylab = expression(T[h](V)))
dev.off()
##################################################

pdf("densityplots.pdf", width = 7, height = 5)
opar <- par(pty = "s", 
            mar = c(5.1, 4.1, 4.1, 2.1) - 2, 
            mgp=c(1.4,0.5,0))
lay <- matrix(1:6, ncol = 3, byrow = TRUE) # layout matrix
layout(lay)
u <- seq(from=0, to = 1, length=40)
v <- u
Tu <- udpbex(u, expansion1, 1)
Tv <- udpbex(v, expansion1, 2)
cop <- outer(Tu, Tv, FUN = function(x,y){dCopula(cbind(x,y), copmod1@copula)})
persp(u, v, cop, zlab = "c(u,v)")
Tu <- udpbex(u, expansion2, 1)
Tv <- udpbex(v, expansion2, 2)
cop <- outer(Tu, Tv, FUN = function(x,y){dCopula(cbind(x,y), copmod2@copula)})
persp(u, v, cop, zlab = "c(u,v)")
Tu <- udpbex(u, expansion3, 1)
Tv <- udpbex(v, expansion3, 2)
cop <- outer(Tu, Tv, FUN = function(x,y){dCopula(cbind(x,y), copmod3@copula)})
persp(u, v, cop, zlab = "c(u,v)")
cop <- outer(u, v, FUN = function(x,y){copdens(cbind(x,y), model = "Gumbel", theta = parresults1$par.ests, case=1)})
persp(u, v, cop, zlab = "c(u,v)")
cop <- outer(u, v, FUN = function(x,y){copdens(cbind(x,y), model = "Gumbel", theta = parresults2$par.ests, case=2)})
persp(u, v, cop, zlab = "c(u,v)")
cop <- outer(u, v, FUN = function(x,y){copdens(cbind(x,y), model = "Frank", theta = parresults3$par.ests, case=3)})
persp(u, v, cop, zlab = "c(u,v)")
par(opar)
dev.off()