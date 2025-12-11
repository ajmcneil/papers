# Code independent of basiscor library

# Case of cosine basis

cos_solve <- function(x, n = 1){
  theta <- acos(x)
  sign <- ((-1)^(1:n)) * (-1)
  top <- floor(n/2)
  add <- c(0, rep(1:top, each = 2))[1:n]
  sol <- (theta*sign + add*2*pi)/(n*pi)
  cbind(sol, sign, add)
}
cos_solve(0.5, 7)

cosplot <- function(j, k, mincop = FALSE, np =500){
  uvals <- seq(from = 0, to = 1, length = np)
  for (i in 1:length(uvals)){
    rhs <- (-1)^(j + k) * cos(j * pi * uvals[i])
    if (mincop)
      rhs <- rhs * -1
    v <- cos_solve(rhs, k)[,1]
    if (i == 1)
      plot(rep(uvals[1], length(v)), v, 
       xlim = c(0,1), ylim = c(0,1), axes=FALSE, xlab="", ylab="", cex=0.3)
    else
      points(rep(uvals[i], length(v)), v, cex = 0.3)
  }
}


par(mfrow = c(5,5))
for (j in 1:5)
  for (k in 1:5)
    cosplot(j,k)
for (j in 1:5)
  for (k in 1:5)
    cosplot(j,k, TRUE)
par(mfrow=c(1,1))

jmax <- 6
pdf("cosmaxplots.pdf")
par(mfrow = c(jmax,jmax),     
    oma = c(1, 1, 1, 1), 
    mar = c(1,1, 0, 0), 
    mgp = c(0, 0, 0),    
    xpd = NA)  
for (j in 1:jmax)
  for (k in 1:jmax){
    cosplot(j, k)
    axis(side=1, labels=FALSE, tick = FALSE)
    axis(side=2, labels=FALSE, tick = FALSE)
    box()
  }
dev.off()

jmax <- 6
pdf("cosminplots.pdf")
par(mfrow = c(jmax,jmax),     
    oma = c(1, 1, 1, 1), 
    mar = c(1,1, 0, 0), 
    mgp = c(0, 0, 0),    
    xpd = NA)  
for (j in 1:jmax)
  for (k in 1:jmax){
    cosplot(j, k, TRUE)
    axis(side=1, labels=FALSE, tick = FALSE)
    axis(side=2, labels=FALSE, tick = FALSE)
    box()
  }
dev.off()

# Case of Legendre basis

plotextremalLegendre <- function(j,
                                 k,
                                 type = "max",                      
                                 n = 5000,
                                 nsim = 10000,
                                 gridlines = TRUE,
                                 pointsize = 0.3,
                                 ...) {
  uvalues <- seq(from = 0,
                 to = 1,
                 length = n)
  eps <- 1e-06
  polyj <- orthopolynom::slegendre.polynomials(j + 1)[[j + 1]]
  polyk <- orthopolynom::slegendre.polynomials(k + 1)[[k + 1]]
  Vr <- matrix(nrow = n, ncol = k)
  Vi <- Vr
  dosims <- FALSE
  if ((k != j) | (type == "min"))
  {
    dosims <- TRUE
    Usim <- runif(nsim)
    PUj <- predict(polyj, Usim)
    PUk <- predict(polyk, Usim)
    Fj <- ecdf(PUj)
  }
  lhs <- predict(polyj, uvalues)
  if (dosims) {
    ppts <- Fj(lhs)
    if (type == "min")
      ppts <- 1- ppts
    lhs <- quantile(PUk, ppts)
  }
  cfs <- coef(polyk)
  for (i in 1:length(uvalues)) {
    newcfs <- cfs
    newcfs[1] <- cfs[1] - lhs[i]
    roots <- polyroot(newcfs)
    Vr[i, ] <- Re(roots)
    Vi[i, ] <- Im(roots)
  }
  select <- abs(Vi[, 1]) <= eps
  U <- uvalues[select]
  V <- Vr[select, 1]
  plot(
    U,
    V,
    xlim = c(0, 1),
    ylim = c(0, 1), 
    cex = pointsize,
    ...
  )
  if (k > 1) {
    for (i in 2:k) {
      select <- abs(Vi[, i]) <= eps
      points(uvalues[select], Vr[select, i], cex= pointsize)
    }
    if (gridlines) {
      tpsk <- polyroot(coef(deriv(polyk, "x")))
      for (i in 1:length(tpsk))
        abline(h = Re(tpsk[i]),
               col = "red",
               lty = 2)
    }
  }
  if ((j > 1) & gridlines) {
    tpsj <- polyroot(coef(deriv(polyj, "x")))
    for (i in 1:length(tpsj))
      abline(v = Re(tpsj[i]),
             col = "red",
             lty = 2)
  }
}

par(pty = "s", mfrow=c(1,1),
    mgp = c(2, 1, 0))
plotextremalLegendre(1,1)
plotextremalLegendre(2,2)
plotextremalLegendre(4, 4, "min")
plotextremalLegendre(4,4)
plotextremalLegendre(3, 3)
plotextremalLegendre(2,5, "min")
plotextremalLegendre(4,3, "min")
plotextremalLegendre(5,7, "min") 

jmax <- 6

pdf("maxplots.pdf")
par(mfrow = c(jmax,jmax),     
    oma = c(1, 1, 1, 1), 
    mar = c(1,1, 0, 0), 
    mgp = c(0, 0, 0),    
    xpd = NA)  
for (j in 1:jmax)
  for (k in 1:jmax){
    plotextremalLegendre(j, k, n=10000, gridlines=FALSE, axes=FALSE, xlab="", ylab="")
    axis(side=1, labels=FALSE, tick = FALSE)
    axis(side=2, labels=FALSE, tick = FALSE)
    box()
  }
dev.off()


pdf("minplots.pdf")
par(mfrow = c(jmax,jmax),     
    oma = c(1, 1, 1, 1), 
    mar = c(1,1, 0, 0), 
    mgp = c(0, 0, 0),    
    xpd = NA)  
for (j in 1:jmax)
  for (k in 1:jmax){
    plotextremalLegendre(j, k, "min", n=10000, gridlines=FALSE, axes=FALSE, xlab="", ylab="")
    axis(side=1, labels=FALSE, tick = FALSE)
    axis(side=2, labels=FALSE, tick = FALSE)
    box()
  }
dev.off()

# Focus on (3,3) and (4,4) cases of Legendre

n <- 10000
U <- seq(from=0,to =1, length=n)
eps <- 1e-6


Vr <- matrix(nrow=n, ncol = 3)
Vi <- Vr
for (i in 1:length(U)){
  lhs <- 20*U[i]^3 - 30*U[i]^2 + 12*U[i]
  roots <- polyroot(c(-lhs, 12, -30, 20))
  Vr[i,] <- Re(roots)
  Vi[i,] <- Im(roots)
}
pdf("P3_max.pdf")
par(pty = "s", mfrow=c(1,1), mgp =c(2,1,0))
select1 <- abs(Vi[,1]) <= eps
plot(U[select1], Vr[select1,1], xlim= c(0,1), ylim=c(0,1),xlab="U", ylab="V")
select2 <- abs(Vi[,2]) <= eps
points(U[select2], Vr[select2,2])
select3 <- abs(Vi[,3]) <= eps
points(U[select3], Vr[select3,3])


tps <- polyroot(c(12,-60,60))
abline(h=Re(tps[1]),col="red", lty=2)
abline(h=Re(tps[2]),col="red", lty=2)
abline(v=Re(tps[1]),col="red", lty=2)
abline(v=Re(tps[2]),col="red", lty=2)
dev.off()

Vr <- matrix(nrow=n, ncol = 4)
Vi <- Vr
for (i in 1:length(U)){
  lhs <- 70*U[i]^4 - 140*U[i]^3 + 90*U[i]^2 -20 *U[i]
  roots <- sort(polyroot(c(-lhs, -20, 90, -140, 70)))
  Vr[i,] <- Re(roots)
  Vi[i,] <- Im(roots)
}
pdf("P4_max.pdf")
par(pty = "s", mfrow=c(1,1), mgp =c(2,1,0))
select1 <- abs(Vi[,1]) <= eps
plot(U[select1],Vr[select1,1],xlim= c(0,1),ylim=c(0,1),xlab="U",ylab="V")
select2 <- abs(Vi[,2]) <= eps
points(U[select2],Vr[select2,2])
select3 <- abs(Vi[,3]) <= eps
points(U[select3],Vr[select3,3])
select4 <- abs(Vi[,4]) <= eps
points(U[select4],Vr[select4,4])

tps <- polyroot(c(-20,180,-420, 280))
abline(h=Re(tps[1]),col="red", lty=2)
abline(h=Re(tps[2]),col="red", lty=2)
abline(h=Re(tps[3]),col="red", lty=2)
abline(v=Re(tps[1]),col="red", lty=2)
abline(v=Re(tps[2]),col="red", lty=2)
abline(v=Re(tps[3]),col="red", lty=2)
dev.off()

pdf("P4_max_vt.pdf")
par(pty = "s", mfrow=c(1,1), mgp =c(2,1,0))
plot(abs(2*U[select1]-1),abs(2*Vr[select1,1]-1),xlim= c(0,1),ylim=c(0,1),xlab="U",ylab="V")
points(abs(2*U[select2]-1),abs(2*Vr[select2,2]-1))
dev.off()

pdf("P4_max2.pdf")
par(pty = "s", mfrow=c(1,1), mgp =c(2,1,0))
theta <- runif(n)*2*pi
r <- sqrt(3/14)
u <- cos(theta)*r + 0.5
v <- sin(theta)*r + 0.5
plot(u, v, ylim = c(0,1), xlim = c(0,1),xlab="U",ylab="V")
u1 <- 0.5 - sqrt(r^2/2)
u2 <- 0.5 + sqrt(r^2/2)
lines(c(u1,u2), c(u2,u1), lwd =10)
lines(c(0,0.5-r), c(0, 0.5-r),lwd = 10)
lines(c(r+0.5,1), c(r+0.5, 1),lwd = 10)
tps <- polyroot(c(-20,180,-420, 280))
abline(h=Re(tps[1]),col="red", lty=2)
abline(h=Re(tps[2]),col="red", lty=2)
abline(h=Re(tps[3]),col="red", lty=2)
abline(v=Re(tps[1]),col="red", lty=2)
abline(v=Re(tps[2]),col="red", lty=2)
abline(v=Re(tps[3]),col="red", lty=2)
dev.off()
