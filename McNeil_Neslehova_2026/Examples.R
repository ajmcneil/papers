library(udp)
library(rvinecopulib)

#####################################################################
# Motivating example (Figure 1)

set.seed(13)
n <- 1000
cop <- bicop_dist("gumbel", parameters = 2.5)
V <- rbicop(n, cop)
V1 <- V[,1]
V2 <- V[,2]
U1 <- cbind(udpsi(vsymmetric(), V1), udpsi(vsymmetric(), V2))
U2 <- cbind(udpsi(udpcosine(3), V1), udpsi(udpid(), V2))
U3 <- cbind(udpsi(udpcosine(3), V1), udpsi(udpcosine(3), V2))
X1 <- apply(U1,2,qnorm)
X2 <- apply(U2,2,qnorm)
X3 <- apply(U3,2,qnorm)
xr <- range(X1,X2,X3)

pdf("motivationA.pdf",width=8, height=2)
par(mfrow = c(1,4),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")
plot(X1[,1], X1[,2], xlim =xr, ylim = xr, xlab=expression(X), ylab =expression(Y))
plot(X2[,1], X2[,2], xlim =xr, ylim = xr, xlab=expression(X), ylab =expression(Y))
plot(X3[,1], X3[,2], xlim =xr, ylim = xr, xlab=expression(X), ylab =expression(Y))
plot(qchisq(V1,1),qchisq(V2,1), xlab = expression(g(X)),ylab=expression(h(Y)))
dev.off()
pdf("motivationB.pdf",width=8, height=2)
par(mfrow = c(1,4),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")
plot(U1[,1], U1[,2], xlim = c(0,1), ylim = c(0,1), xlab=expression(U[1]), ylab =expression(U[2]))
plot(U2[,1], U2[,2], xlim = c(0,1), ylim = c(0,1), xlab=expression(U[1]), ylab =expression(U[2]))
plot(U3[,1], U3[,2], xlim = c(0,1), ylim = c(0,1), xlab=expression(U[1]), ylab =expression(U[2]))
plot(V1, V2, xlim = c(0,1), ylim= c(0,1), xlab = expression(V[1]),ylab=expression(V[2]))
dev.off()


x <- seq(from = -3, to = 3, length=400)
g <- function(x){qchisq(udptrans(vsymmetric(), pnorm(x)), 1)}
plot(x, g(x), type="l")

pdf("supplement.pdf",width=6.5, height=3)
par(mfrow = c(1,2),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")
plot(udpcosine(3))
h <- function(x){qchisq(udptrans(udpcosine(3),pnorm(x)), 1)}
plot(x, h(x), type="l", ylab= "g(x)")
dev.off()

#############################################################
# Regular udp functions (Figure 2)

pdf("TplotsLegendre.pdf",width=8, height=2)
par(mfrow = c(1,5),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")
for (j in 2:6)
  plot(udplegendre(j), embellish = "bw")
par(mfrow = c(1,1))
dev.off()

##################################################################
# Bivariate examples from Gauss copula and v-transforms (Figure 3)

n <- 100
u1 <- seq(from = 0.01, to = 0.99, length =n)
u2 <- u1
vt1 <- vsymmetric()
vt2 <- vsymmetric()
copV <- bicop_dist("gauss", parameters = 0.85)
mod <- bsicopula(copV, vt1,vt2)
z <- outer(u1, u2, FUN = dbsicopula, object = mod)
sum(z)/(n^2)
pdf("Gauss1.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
#################################################################################

wtmodel <- randsdvine(bicop_dist("gauss", parameters = 0.8),
                      bicop_dist("gauss", parameters = 0),
                      bicop_dist("gauss", parameters = 0))
mod <- bsicopula(copV, vt1, vt2, wtmodel)
z <- outer(u1, u2, FUN = dbsicopula, object = mod)
sum(z)/(n^2)
pdf("Gauss2.pdf")
contour(u1, u2, z, nlevel = 100)
dev.off()
#####################################################################

wtmodel <- randmixture(cop1 = bicop_dist("gauss", parameters = 1),
                       cop2 = bicop_dist("gauss", parameters = -1),
                       function(v1, v2){pmax(v1, v2) > 0.6})
mod <- bsicopula(copV, vt1, vt2, wtmodel)
z <- outer(u1, u2, FUN = dbsicopula, object = mod)
sum(z)/(n^2)
pdf("Gauss3.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
######################################################################

wtmodel <- randsdvine(bicop_dist("gauss", parameters = 0.8),
                      bicop_dist("gauss", parameters = 0.7),
                      bicop_dist("gauss", parameters = 0.1))
mod <- bsicopula(copV, vt1, vt2, wtmodel)
z <- outer(u1, u2, FUN = dbsicopula, object = mod)
sum(z)/(n^2)
z[z < 0.00001] <- 0
pdf("Gauss4.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
