library(udp)
library(rvinecopulib)
library(qrmdata)

data("SP500_const", package = "qrmdata")
stock <- SP500_const["2006-01-01/2015-12-31", "JPM"]
plot(stock)
X <- (diff(log(stock))[-1])*100 # log-returns
X <- X[X!=0] # 0 log-returns from market closed days
pdf("tsplot.pdf",width=6, height=4)
plot(X)
dev.off()

U <- as.numeric(pseudo_obs(X))
data <- cbind(U[-length(U)], U[-1])
pdf("tsdataU.pdf",width=4, height=4)
par(mar = c(2, 3, 0, 0),
    mgp = c(1.2, 0.5, 0))
plot(data, xlab = expression(U[t]), ylab = expression(U[t+1]))
dev.off()

mod1 <- bsicopula(bicop_dist(family = "joe", parameters = 1.5),
                  udp1 = vlinear(0.5),
                  udp2 = vlinear(0.5))
(fit1 <- fitbsicopula(data, mod1, se = "hessian"))

randmodA <- randsdvine(bicop_dist(family = "gauss", parameters = 0))
mod2 <- bsicopula(bicop_dist(family = "joe", parameters = 1.5),
                  udp1 = vlinear(0.5),
                  udp2 = vlinear(0.5),
                  randomizermod = randmodA)
(fit2 <- fitbsicopula(data, mod2, se = "hessian"))
(p1 <- 1 - pchisq(2*(logLik(fit2) - logLik(fit1)), 1))

randmodB <- randsdvine(bicop_dist(family = "gauss", parameters = -0.11),
                       bicop_dist(family = "gauss", parameters = 0),
                       bicop_dist(family = "gauss", parameters = 0))
mod3 <- bsicopula(bicop_dist(family = "joe", parameters = 1.5),
                  udp1 = vlinear(0.5),
                  udp2 = vlinear(0.5),
                  randomizermod = randmodB)
(fit3 <- fitbsicopula(data, mod3, se = "hessian"))
(p2 <- 1 - pchisq(2*(logLik(fit3) - logLik(fit2)), 2))

pdf("tsdataV.pdf", width=4, height=4)
par(mar = c(2, 3, 0, 0),
    mgp = c(1.2, 0.5, 0))
V1 <- udptrans(fit1@bsicopula@udp1, data[,1])
V2 <- udptrans(fit1@bsicopula@udp2, data[,2])
plot(cbind(V1, V2), xlab = expression(V[t]), ylab = expression(V[t+1]))
dev.off()

