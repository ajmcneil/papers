library(udp)

set.seed(13)
n <- 1000
cop <- copula::gumbelCopula(2.5)
V <- copula::rCopula(n,cop)
V1 <- V[,1]
V2 <- V[,2]
plot(V1,V2)
U1 <- cbind(udpsi(vsymmetric(), V1), udpsi(vsymmetric(), V2))
U2 <- cbind(udpsi(udpcosine(3), V1), udpsi(shuffle(1), V2))
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
