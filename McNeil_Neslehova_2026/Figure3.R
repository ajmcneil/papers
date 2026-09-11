library(udp)
library(rvinecopulib)
source("functions.R")




############################################################################

n <- 100
u1 <- seq(from = 0.01, to = 0.99, length =n)
u2 <- u1
vt1 <- vlinear(0.5)
vt2 <- vlinear(0.5)
copV <- rvinecopulib::bicop_dist("gauss", parameters = 0.85)
z <- outer(u1, u2, FUN = bsicopula, udp1 = vt1, udp2 = vt2,
           copV = copV)
sum(z)/(n^2)
pdf("Gauss1.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
#################################################################################

wtmodel <- list(type = "dvine",
                copZ1V2 = bicop_dist("gauss", parameters = 0),
                copV1Z2 = bicop_dist("gauss", parameters = 0),
                copZV = bicop_dist("gauss", parameters = 0.8))
z <- outer(u1, u2, FUN = bsicopula, udp1 = vt1, udp2 = vt2,
           copV = copV, wtmodel)
sum(z)/(n^2)
pdf("Gauss2.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
#####################################################################

wtmodel <- list(type = "condcopula",
                switch = function(v1, v2, k){
                  pmax(v1, v2) > k
                },
                k = 0.6,
                cop1 = bicop_dist("gauss", parameters = 1),
                cop2 = bicop_dist("gauss", parameters = -1))
z <- outer(u1, u2, FUN = bsicopula, udp1 = vt1, udp2 = vt2,
           copV = copV, wtmodel)
sum(z)/(n^2)
apply(z,1,sum)/n
apply(z,2,sum)/n
pdf("Gauss3.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
######################################################################

wtmodel <- list(type = "dvine",
                copZ1V2 = bicop_dist("gauss", parameters = 0.7),
                copV1Z2 = bicop_dist("gauss", parameters = 0.1),
                copZV = bicop_dist("gauss", parameters = 0.8))
z <- outer(u1, u2, FUN = bsicopula, udp1 = vt1, udp2 = vt2,
           copV = copV, wtmodel = wtmodel)
sum(z)/(n^2)
pdf("Gauss4.pdf")
contour(u1, u2, z,nlevel = 100)
dev.off()
