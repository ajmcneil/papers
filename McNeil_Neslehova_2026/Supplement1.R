library(udp)
library(rvinecopulib)
source("functions.R")

copV <- bicop_dist("gauss", parameters = 0.85)
vt1 <- vlinear(0.5)
vt2 <- vlinear(0.5)
wtmodel <- list(type = "dvine", 
                copZ1V2 = bicop_dist("gauss", parameters = 0.7),
                copV1Z2 = bicop_dist("gauss", parameters = 0.1),
                copZV = bicop_dist("gauss", parameters = 0.8))

FZgivVmarg <- function(z, v1, v2, copV, wtmodel, marg){
  if (marg == 1){
    arg1 <- hbicop(cbind(v1, v2), cond_var = 1, copV)
    return(hbicop(cbind(z, arg1), cond_var = 2, wtmodel$copZ1V2))
  }
  if (marg == 2){
    arg2 <- hbicop(cbind(v1, v2), cond_var = 2, copV)
    return(hbicop(cbind(z, arg2), cond_var = 2, wtmodel$copV1Z2))
  }
}
integrand <- function(v, u, copV, wtmodel, marg){
  if (marg == 1)
    return(FZgivVmarg(0.5, abs(2*u-1), v, copV, wtmodel, marg))
  if (marg == 2)
    return(FZgivVmarg(0.5, v, abs(2*u-1), copV, wtmodel, marg))
}
# check integrand
u <- 0.1
v <- seq(from = 0, to = 1, length =100)
Fv <- integrand(v, u, copV, wtmodel, 1)
plot(v,Fv)

# pdf("margdens.pdf")
# par(mfrow=c(1,1),cex=1.5)
# uvals <- seq(from = 0, to = 1, length = 200)
# omega1 <- rep(NA,length(uvals))
# for (i in 1:length(uvals)){
#   if (uvals[i] < 0.5)
#    omega1[i] <- 2*integrate(integrand, lower = 0, upper = 1, u= uvals[i], copV = copV, wtmodel = wtmodel, marg = 1)$value
#   else
#     omega1[i] <- 2*(1-integrate(integrand, lower = 0, upper = 1, u= uvals[i], copV = copV, wtmodel = wtmodel, marg = 1)$value) 
# }
# plot(uvals, omega1, type="l", xlab = "u", ylab = expression(omega[1](u)))
# dev.off()

pdf("margdens.pdf")
par(mfrow=c(1,1),cex=1.5)
uvals1 <- seq(from = 0, to = 0.499, length = 100)
omega11 <- rep(NA,length(uvals1))
uvals2 <- seq(from = 0.501, to = 1, length = 100)
omega12 <- rep(NA,length(uvals2))
for (i in 1:length(uvals1)){
    omega11[i] <- 2*integrate(integrand, lower = 0, upper = 1, u= uvals1[i], copV = copV, wtmodel = wtmodel, marg = 1)$value
    omega12[i] <- 2*(1-integrate(integrand, lower = 0, upper = 1, u= uvals2[i], copV = copV, wtmodel = wtmodel, marg = 1)$value) 
}
plot(NULL, xlim = c(0,1), ylim = c(0,2),xlab = "u", ylab = expression(omega[1](u)))
lines(uvals1, omega11)
lines(uvals2, omega12)
dev.off()

n <- 100
u1 <- seq(from = 0.01, to = 0.99, length =n)
u2 <- u1
z <- outer(u1, u2, FUN = bsicopula, udp1 = vt1, udp2 = vt2, 
           copV = copV, wtmodel = wtmodel, wtonly = TRUE)
sum(z)/(n^2)
contour(u1, u2, z,nlevel = 100)
persp(u1,u2,z, zlab = "w(u1, u2)", theta =240, phi = 30) 
persp(u1,u2,z, zlab = "w(u1, u2)", theta =140, phi = 30) 
plot(u1,apply(z,1,sum)/n)
plot(u2,apply(z,2,sum)/n)
##############################################################


