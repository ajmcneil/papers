library(udp)
library(rvinecopulib)
library(lattice)
library(gridExtra)
source("functions.R")

copV <- bicop_dist("gauss", parameters = 0.85)
vt1 <- vlinear(0.5)
vt2 <- vlinear(0.5)
wtmodel <- list(type = "dvine", 
                copZ1V2 = bicop_dist("gauss", parameters = 0.7),
                copV1Z2 = bicop_dist("gauss", parameters = 0.1),
                copZV = bicop_dist("gauss", parameters = 0.8))


multinomialp <- function(v1, v2, delta1, delta2, copV, wtmodel = NULL, quadrant = 1){
  switch(quadrant,
         FZgivV(delta1, delta2, v1, v2, copV, wtmodel),
         FZgivV(delta1, 1, v1, v2, copV, wtmodel) - FZgivV(delta1, delta2, v1, v2, copV, wtmodel),
         FZgivV(1, delta2, v1, v2, copV, wtmodel) - FZgivV(delta1, delta2, v1, v2, copV, wtmodel),
         1 - FZgivV(1, delta2, v1, v2, copV, wtmodel) - FZgivV(delta1, 1, v1, v2, copV, wtmodel) +
           FZgivV(delta1, delta2, v1, v2, copV, wtmodel)
  )
}



#####################################################################

integrand <- function(v2, v1, delta1, delta2, copV, wtmodel, quadrant){
  multinomialp(v1, v2, delta1, delta2, copV, wtmodel, quadrant)*
    dbicop(cbind(v1, v2), family = copV)
}

v1 <- seq(from = 0.001, to = 0.999, length = 100)
pr <- matrix(nrow = length(v1), ncol =4, NA)
for (i in 1: length(v1)){
pr[i,1] <- integrate(integrand,0, 1, v1 = v1[i], delta1 = 0.5, delta2 = 0.5, 
                copV = copV, wtmodel = wtmodel, quadrant =1)$value
pr[i,2] <- integrate(integrand,0, 1, v1 = v1[i], delta1 = 0.5, delta2 = 0.5, 
                copV = copV, wtmodel = wtmodel, quadrant =2)$value
pr[i,3] <- integrate(integrand,0, 1, v1 = v1[i], delta1 = 0.5, delta2 = 0.5, 
                copV = copV, wtmodel = wtmodel, quadrant =3)$value
pr[i,4] <- integrate(integrand,0, 1, v1 = v1[i], delta1 = 0.5, delta2 = 0.5, 
                copV = copV, wtmodel = wtmodel, quadrant =4)$value}
pr[,1]+pr[,2]
pr[,3]+pr[,4]
pr[,2]+pr[,4]


pdf("quadplot1.pdf", width = 5, height = 5)
par(mfrow=c(1,1))
plot(v1, pr[,1],
     type = "n",
     xlim = c(0,1),
     ylim = c(0, 1),
     xlab = expression(v),
     ylab = "",
     xaxs = "i", yaxs = "i",
     xaxt = "n")

axis(1)

# Cumulative sums
F1 <- pr[,1]
F2 <- F1 + pr[,2]
F3 <- F2 + pr[,3]
F4 <- F3 + pr[,4]

# Colours
cols <- c("lightblue", "lightgreen", "lightpink", "lightyellow")

polygon(c(v1, rev(v1)),
        c(rep(0, length(v1)), rev(F1)),
        col = cols[1], border = NA)

polygon(c(v1, rev(v1)),
        c(F1, rev(F2)),
        col = cols[2], border = NA)

polygon(c(v1, rev(v1)),
        c(F2, rev(F3)),
        col = cols[3], border = NA)

polygon(c(v1, rev(v1)),
        c(F3, rev(F4)),
        col = cols[4], border = NA)

lines(v1, F1)
lines(v1, F2)
lines(v1, F3)
lines(v1, F4)


# legend(x=0.6, y=1,,
#        legend = c(expression(p[1 * "," * 1](v[1])), 
#                   expression(p[1 * "," * 2](v[1])),
#                   expression(p[2 * "," * 1](v[1])),
#                   expression(p[2 * "," * 2](v[1]))),
#        fill = cols,
#        bty = "n")
legend(x=0.8, y=1,,
       legend = c(expression(Q[1]), 
                  expression(Q[2]),
                  expression(Q[3]),
                  expression(Q[4])),
       fill = cols,
       bty = "n")
dev.off()

######################################################################

integrand <- function(v1, v2, delta1, delta2, copV, wtmodel, quadrant){
  multinomialp(v1, v2, delta1, delta2, copV, wtmodel, quadrant)*
    dbicop(cbind(v1, v2), family = copV)
}

v2 <- seq(from = 0.001, to = 0.999, length = 100)
pr <- matrix(nrow = length(v2), ncol =4, NA)
for (i in 1: length(v2)){
  pr[i,1] <- integrate(integrand,0, 1, v2 = v2[i], delta1 = 0.5, delta2 = 0.5, 
                       copV = copV, wtmodel = wtmodel, quadrant =1)$value
  pr[i,2] <- integrate(integrand,0, 1, v2 = v2[i], delta1 = 0.5, delta2 = 0.5, 
                       copV = copV, wtmodel = wtmodel, quadrant =2)$value
  pr[i,3] <- integrate(integrand,0, 1, v2 = v2[i], delta1 = 0.5, delta2 = 0.5, 
                       copV = copV, wtmodel = wtmodel, quadrant =3)$value
  pr[i,4] <- integrate(integrand,0, 1, v2 = v2[i], delta1 = 0.5, delta2 = 0.5, 
                       copV = copV, wtmodel = wtmodel, quadrant =4)$value}
pr[,1]+pr[,3]
pr[,2]+pr[,4]

pdf("quadplot2.pdf", width = 5, height = 5)
par(mfrow=c(1,1))
plot(v1, pr[,1],
     type = "n",
     xlim = c(0,1),
     ylim = c(0, 1),
     xlab = expression(v),
     ylab = "",
     xaxs = "i", yaxs = "i",
     xaxt = "n")

axis(1)

# Cumulative sums
F1 <- pr[,1]
F2 <- F1 + pr[,3]
F3 <- F2 + pr[,2]
F4 <- F3 + pr[,4]

# Colours
cols <- c("lightblue", "lightpink","lightgreen", "lightyellow")

polygon(c(v2, rev(v2)),
        c(rep(0, length(v2)), rev(F1)),
        col = cols[1], border = NA)

polygon(c(v2, rev(v2)),
        c(F1, rev(F2)),
        col = cols[2], border = NA)

polygon(c(v2, rev(v2)),
        c(F2, rev(F3)),
        col = cols[3], border = NA)

polygon(c(v2, rev(v2)),
        c(F3, rev(F4)),
        col = cols[4], border = NA)

lines(v2, F1)
lines(v2, F2)
lines(v2, F3)
lines(v2, F4)

# legend(x=0.6, y=1,,
#        legend = c(expression(p[1 * "," * 1](v[2])), 
#                   expression(p[2 * "," * 1](v[2])),
#                   expression(p[1 * "," * 2](v[2])),
#                   expression(p[2 * "," * 2](v[2]))),
#        fill = cols,
#        bty = "n")
legend(x=0.8, y=1,,
       legend = c(expression(Q[1]), 
                  expression(Q[3]),
                  expression(Q[2]),
                  expression(Q[4])),
       fill = cols,
       bty = "n")
dev.off()
########################################################

n <- 100
v <- seq(from = 0, to = 1, length =n)
layout(matrix(1:4, 2, 2, byrow = TRUE))
z2 <- outer(v, v, multinomialp, delta1 =0.5, delta2 = 0.5, copV = copV, wtmodel = wtmodel,
            quadrant = 2)
z4 <- outer(v, v, multinomialp, delta1 =0.5, delta2 = 0.5, copV = copV, wtmodel = wtmodel,
            quadrant = 4)
z1 <- outer(v, v, multinomialp, delta1 =0.5, delta2 = 0.5, copV = copV, wtmodel = wtmodel,
            quadrant = 1)
z3 <- outer(v, v, multinomialp, delta1 =0.5, delta2 = 0.5, copV = copV, wtmodel = wtmodel,
            quadrant = 3)

summary(as.vector(z2 + z4 + z1 + z3))
z <- list(z1, z2, z3, z4)

# Common colour scale
zlim <- range(unlist(z), na.rm = TRUE)
at <- seq(zlim[1], zlim[2], length.out = 11)

ps <- list(
  top.padding = 0,
  bottom.padding = 0,
  left.padding = 0,
  right.padding = 0,
  axis.padding = 0
)

xr <- c(0, 1)

p1 <- levelplot(z1, 
                row.values = v,
                column.values = v,
                xlim = xr,
                ylim = xr,
                xlab = expression(v[1]),
                ylab = expression(v[2]),
           #     main = expression(p[1 * "," * 1](v[1],v[2])),
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
         #       main = expression(p[1 * "," * 2](v[1],v[2])),
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
       #         main = expression(p[2 * "," * 1](v[1],v[2])),
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
      #          main = expression(p[2 * "," * 2](v[1],v[2])),
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