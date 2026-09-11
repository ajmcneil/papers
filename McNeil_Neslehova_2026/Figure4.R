library(udp)
library(rvinecopulib)
library(qrmdata)
source("functions.R")

data("SP500_const", package = "qrmdata")
stock <- SP500_const["2006-01-01/2015-12-31","JPM"]
plot(stock)
X <- (diff(log(stock))[-1])*100 # log-returns
X <- X[X!=0] # 0 log-returns from market closed days

length(X)
U <- rank(X, na.last= "keep", ties.method = "random")/(length(X) + 1)
data <- cbind(U[-length(U)], U[-1])
(fit <- fit_bsicopula(data))
V1 <- udptrans(vlinear(fit$stage1[1]), data[,1])
V2 <- udptrans(vlinear(fit$stage1[2]), data[,2])

par(mfrow = c(1,1),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")

pdf("tsplot.pdf",width=6, height=4)
plot(X)
dev.off()

pdf("tsdataU.pdf",width=4, height=4)
par(mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0))
plot(data, xlab = expression(U[t]), ylab = expression(U[t+1]))
dev.off()

pdf("tsdataV.pdf", width=4, height=4)
par(mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0))
plot(cbind(V1, V2), xlab = expression(V[t]), ylab = expression(V[t+1]))
dev.off()

