library(udp)
library(rvinecopulib)
library(lattice)
library(gridExtra)

# set-up

copV <- bicop_dist("gauss", parameters = 0.85)
vt1 <- vlinear(0.5)
vt2 <- vlinear(0.5)
wtmodel <- randsdvine(bicop_dist("gauss", parameters = 0.8),
                      bicop_dist("gauss", parameters = 0.7),
                      bicop_dist("gauss", parameters = 0.1))
bsimod <- bsicopula(copV, vt1, vt2, wtmodel)

#########################################################################
### Weight function (a density on [0,1]^2)

n <- 100
u1 <- (seq_len(n) - 0.5) / n
u2 <- u1
num <- outer(u1, u2, FUN = dbsicopula, object=bsimod)
den <- outer(u1, u2, FUN = dbsicopula, object=bsicopula(copV, vt1, vt2))
wt <- num/den
persp(u1, u2, wt, zlab = "w(u1, u2)", theta =240, phi = 30)
persp(u1, u2, wt, zlab = "w(u1, u2)", theta =140, phi = 30)

#################################################################
# Marginal distribution of weight function not uniform
# Weight function not a copula

n <- 100
u1 <- seq(from = 0.001, to = 0.4999, length = n)
u2 <- seq(from = 0.5001, to = 0.999, length = n)
u <- c(u1, u2)
v <- (seq_len(2*n) - 0.5) / (2*n)

pc <- pcells(bsimod, v1 = abs(2*u - 1), v2 = v, grid = TRUE)
pcint <- apply(pc,c(1, 2, 3), mean) # numerical integral wrt v
wm <- 2 * apply(pcint, c(1, 3), sum) # sums over cell2 (Z2)
wm1 <- wm[1, ][u < 0.5] # (Z1 < 0.5)
wm2 <- wm[2, ][u > 0.5] # (Z1 > 0.5)

pdf("margdens.pdf")
par(mfrow=c(1,1),cex=1.5)
plot(NULL, xlim = c(0,1), ylim = c(0,2), xlab = "u", ylab = expression(omega[1](u)))
lines(u1, wm1)
lines(u2, wm2)
dev.off()

##############################################################
#### Multinomial probability function

n <- 100
v <- (seq_len(n) - 0.5) / n
G <- pcells(bsimod, v, v, grid = TRUE)
z1 <- G[1,1,,]
z2 <- G[2,1,,]
z3 <- G[1,2,,]
z4 <- G[2,2,,]
z <- list(z1, z2, z3, z4)

zlim <- range(unlist(z), na.rm = TRUE)
at <- seq(zlim[1], zlim[2], length.out = 11)
ps <- list(top.padding = 0, bottom.padding = 0, left.padding = 0, right.padding = 0, axis.padding = 0)
xr <- c(0, 1)

p1 <- levelplot(z1,
                row.values = v,
                column.values = v,
                xlim = xr,
                ylim = xr,
                xlab = expression(v[1]),
                ylab = expression(v[2]),
                main = expression(Q[1]),
                at = at,
                scales = list(
                  x = list(at = seq(0, 1, 0.2)),
                  y = list(at = seq(0, 1, 0.2))
                ),
                par.settings = ps)
p2 <- levelplot(z2,
                row.values = v,
                column.values = v,
                xlim = xr,
                ylim = xr,
                xlab = NULL,
                ylab = expression(v[2]),
                main = expression(Q[2]),
                at = at,
                scales = list(
                  x = list(at = seq(0, 1, 0.2)),
                  y = list(at = seq(0, 1, 0.2))
                ),
                par.settings = ps)
p3 <- levelplot(z3,
                row.values = v,
                column.values = v,
                xlim = xr,
                ylim = xr,
                xlab = expression(v[1]),
                ylab = NULL,
                main = expression(Q[3]),
                at = at,
                scales = list(
                  x = list(at = seq(0, 1, 0.2)),
                  y = list(at = seq(0, 1, 0.2))
                ),
                par.settings = ps)
p4 <- levelplot(z4,
                row.values = v,
                column.values = v,
                xlim = xr,
                ylim = xr,
                xlab = NULL,
                ylab = NULL,
                main = expression(Q[4]),
                at = at,
                scales = list(
                  x = list(at = seq(0, 1, 0.2)),
                  y = list(at = seq(0, 1, 0.2))
                ),
                par.settings = ps)
pdf("levelQ1.pdf", width=5, height=5)
plot(p1)
dev.off()
pdf("levelQ2.pdf", width=5, height=5)
plot(p2)
dev.off()
pdf("levelQ3.pdf", width=5, height=5)
plot(p3)
dev.off()
pdf("levelQ4.pdf", width=5, height=5)
plot(p4)
dev.off()

###################################################################
#### Marginal multinomial probabilities (1)

v <- seq(from = 0.001, to = 0.999, length = 100)
pc <- pcells(bsimod, v1 = v)
F1 <- pc[1,1,]
F2 <- F1 + pc[1,2,]
F3 <- F2 + pc[2,1,]
F4 <- F3 + pc[2,2,]

pdf("quadplot1.pdf", width = 5, height = 5)
par(mfrow=c(1,1))
plot(v, F1,
     type = "n",
     xlim = c(0,1),
     ylim = c(0, 1),
     xlab = expression(v),
     ylab = "",
     xaxs = "i", yaxs = "i",
     xaxt = "n")
axis(1)
cols <- c("lightblue", "lightgreen", "lightpink", "lightyellow")
polygon(c(v, rev(v)),
        c(rep(0, length(v)), rev(F1)),
        col = cols[1], border = NA)
polygon(c(v, rev(v)),
        c(F1, rev(F2)),
        col = cols[2], border = NA)
polygon(c(v, rev(v)),
        c(F2, rev(F3)),
        col = cols[3], border = NA)
polygon(c(v, rev(v)),
        c(F3, rev(F4)),
        col = cols[4], border = NA)
lines(v, F1)
lines(v, F2)
lines(v, F3)
lines(v, F4)
legend(x=0.8, y=1,,
       legend = c(expression(Q[1]),
                  expression(Q[2]),
                  expression(Q[3]),
                  expression(Q[4])),
       fill = cols,
       bty = "n")
dev.off()

###################################################################
#### Marginal multinomial probabilities (2)

v <- seq(from = 0.001, to = 0.999, length = 100)
pc <- pcells(bsimod, v2 = v)
F1 <- pc[1,1,]
F2 <- F1 + pc[2,1,]
F3 <- F2 + pc[1,2,]
F4 <- F3 + pc[2,2,]

pdf("quadplot1.pdf", width = 5, height = 5)
par(mfrow=c(1,1))
plot(v, F1,
     type = "n",
     xlim = c(0,1),
     ylim = c(0, 1),
     xlab = expression(v),
     ylab = "",
     xaxs = "i", yaxs = "i",
     xaxt = "n")
axis(1)
cols <- c("lightblue", "lightpink", "lightgreen", "lightyellow")
polygon(c(v, rev(v)),
        c(rep(0, length(v)), rev(F1)),
        col = cols[1], border = NA)
polygon(c(v, rev(v)),
        c(F1, rev(F2)),
        col = cols[2], border = NA)
polygon(c(v, rev(v)),
        c(F2, rev(F3)),
        col = cols[3], border = NA)
polygon(c(v, rev(v)),
        c(F3, rev(F4)),
        col = cols[4], border = NA)
lines(v, F1)
lines(v, F2)
lines(v, F3)
lines(v, F4)
legend(x=0.8, y=1,,
       legend = c(expression(Q[1]),
                  expression(Q[3]),
                  expression(Q[2]),
                  expression(Q[4])),
       fill = cols,
       bty = "n")
dev.off()

