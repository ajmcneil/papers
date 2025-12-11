library(basiscor)
library(copula)
library(corrplot)

# Basis functions

fourierbasis <- function(u, degree =2){
  if (degree %% 2 == 0)
    output <- sqrt(2)*sin(pi*degree*u)
  else
    output <- sqrt(2)*cos(pi*(degree+1)*u)
  output
}

pdf("basis.pdf", width=7, height=4)
par(mfrow = c(3,5),     # 
    oma = c(1, 1, 1, 1), 
    mar = c(2, 3, 0, 0), 
    mgp = c(0, 1, 0),    
    xpd = NA)  
jrange <- 1:5
u <- seq(from=0.0001, to = 0.9999, length = 100)
for (j in jrange)
  plot(u, basisfunc(u, j, "legendre"), type = "l", ylab = "")
for (j in jrange)
  plot(u, fourierbasis(u, j), type = "l", ylab = "")
for (j in jrange)
  plot(u, basisfunc(u, j, "cosine"), type = "l", ylab = "")
dev.off()

# Chessboard pictures

maxk <- 8
tcop <- matrix(NA,nrow=maxk,ncol=maxk)
sphtcop <- tcop
for (j in 1:maxk)
  for (k in 1:maxk)
    tcop[j,k] <- basiscor(tCopula(0.7, df = 2), j, k)

for (j in 1:maxk)
  for (k in 1:maxk)
    sphtcop[j,k] <- basiscor(tCopula(param = 0, df = 2), j, k)

pdf("radially_symmetric.pdf")
corrplot(tcop, method = "color")
dev.off()
pdf("jointly_symmetric.pdf")
corrplot(sphtcop, method = "color")
dev.off()

# Udp functions

pdf("FplotsLegendre.pdf",width=8, height=4)
par(mfrow = c(2,5),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  
jrange <- 2:6
for (j in jrange){
  Ppoly <- orthopolynom::slegendre.polynomials(j+1)[[j+1]]
  if (j %% 2 != 0)
    lbound <- predict(Ppoly,0)
  else{
    tps <- Re(polyroot(coef(deriv(Ppoly, "x"))))
    lbound <- min(predict(Ppoly, tps)) + 0.002
  }
  x <- seq(from=lbound, to = 1, length=200)
  y <- distLegendre(x, j)$df
  plot(x*sqrt(2*j+1),y, xlab = "x", ylab = expression(F[Lambda[j](U)](x)),type="l")
}
for (j in jrange){
  Ppoly <- orthopolynom::slegendre.polynomials(j+1)[[j+1]]
  if (j %% 2 != 0)
    lbound <- predict(Ppoly,0) 
  else{
    tps <- Re(polyroot(coef(deriv(Ppoly, "x"))))
    lbound <- min(predict(Ppoly, tps)) + 0.002
  }
  x <- seq(from=lbound, to = 1, length=200)
  y <- distLegendre(x, j)$density
  plot(x*sqrt(2*j+1), y/sqrt(2*j+1), xlab = "x", ylab = expression(f[Lambda[j](U)](x)),type="l")
}
par(mfrow = c(1,1))
dev.off()

pdf("TplotsLegendre.pdf",width=8, height=2)
par(mfrow = c(1,5),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpfunc(u, j, "legendre")
  plot(u,v, ylab = expression(T[j]^{Lambda}*(u)), type="l")
}
par(mfrow = c(1,1))
dev.off()


pdf("TplotsCosine.pdf",width=8, height=2)
par(mfrow = c(1,5),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpfunc(u, j, "cosine")
  plot(u,v, ylab = expression(T[j]^{Omega}*(u)), type="l")
}
par(mfrow = c(1,1))
dev.off()

# Check on derivatives

par(mfrow = c(2,5),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpfunc(u, j, "legendre")
  plot(u,v, ylab = expression(T[j]^{Lambda}*(u)), type="l")
}
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpderiv(u, j, "legendre")
  plot(u,v, ylab = expression(T[j]^{Lambda}*(u)), type="l")
}
par(mfrow = c(1,1))

par(mfrow = c(2,5),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpfunc(u, j, "cosine")
  plot(u,v, ylab = expression(T[j]^{Lambda}*(u)), type="l")
}
for (j in jrange){
  u <- seq(from=0, to = 1, length=200)
  v <- udpderiv(u, j, "cosine")
  plot(u,v, ylab = expression(T[j]^{Lambda}*(u)), type="l")
}
par(mfrow = c(1,1))

# Legendre bounds

maxvals <- matrix(NA, nrow= 6, ncol = 6)
diag(maxvals) <- rep(1, 6)
for (j in 1:5)
  for (k in (j+1):6)
    maxvals[j, k] <- extremalLegendre(j,k)
round(maxvals, digits = 3)
#print(xtable::xtable(maxvals, digits = 3), booktabs = TRUE)

minvals <- matrix(NA, nrow= 6, ncol = 6)
for (j in 1:6)
  for (k in j:6)
    minvals[j, k] <- extremalLegendre(j,k, "min")
round(minvals, digits = 3)
#print(xtable::xtable(minvals, digits = 3), booktabs = TRUE)

round(maxvals - minvals, digits = 3)

# Stochastic inversion examples

n <- 1000
set.seed(13)
datastar <- rCopula(n, normalCopula(0.85))
plot(datastar)
j <- 2
k <- 2

pdf("SIpicturee.pdf",width=6, height=2)
par(mfrow = c(1,3),     # 
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 2, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(0, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")  

UV <- rstochinv(datastar, j, k)
plot(UV, xlab = "U", ylab = "V")

Z <- runif(n)
UV <- rstochinv(datastar, j, k, rand = cbind(Z, Z))
plot(UV, xlab = "U", ylab = "V")

Z2 <- rep(NA, n)
for (i in 1:n){
  if (max(datastar[i,1], datastar[i,2]) > 0.6)
    Z2[i] <- Z[i]
  else
    Z2[i] <- 1 - Z[i]
}
UV <- rstochinv(datastar, j, k, rand = cbind(Z, Z2))
plot(UV, xlab = "U", ylab = "V")
dev.off()

# Estimators

n <- 100
set.seed(13)
data <- rCopula(n,claytonCopula(1))

j <- 1
k <- 1
cor(data,method = "spearman")[1,2]
basiscor(data,j,k,method="T1")
basiscor(data,j,k,method="T2")
basiscor(data,j,k,method="T3")
basiscor(data,j,k,method="T4")
basiscor(data,j,k,method="T5")
basiscor(data,j,k,method="T6")

j <- 1
k <- 2
basiscor(data,j,k,method="T1")
basiscor(data,j,k,method="T2")
basiscor(data,j,k,method="T3")
basiscor(data,j,k,method="T4")
basiscor(data,j,k,method="T5")
basiscor(data,j,k,method="T6")

j <- 2
k <- 2
cor((apply(data,2,rank)/(n+1) - 0.5)^2)
cor((2*apply(data,2,rank)-n-1)^2)
basiscor(data,j,k,method="T1")
basiscor(data,j,k,method="T2")
basiscor(data,j,k,method="T3")
basiscor(data,j,k,method="T4")
basiscor(data,j,k,method="T5")
basiscor(data,j,k,method="T6")