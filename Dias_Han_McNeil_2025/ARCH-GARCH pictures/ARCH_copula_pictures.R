require(skewt)
source("Section2_functions.R")

# ARCH(1) parameters
ARCHpars <- c(0.4, 0.6, 0)
# ARCH(1) leverage parameters
ARCHlevpars <- c(0.4, 0.3, 0)
levpar <- 0.4
# AR(1) -ARCH(1)
phi <- 0.3
# Parameters of innovation distributions
nuA <- 4
nuB <- 2.5
lambda <- 0.8



#####################################################
# ARCH(1) with Gaussian innovations
#####################################################

set.seed(123)
X <- simgarch(dist = "norm", n=10000, pars = ARCHpars, copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_norm <- marg_approx(c(-5,5), innov = "norm", pars = ARCHpars)
marg_norm$ev_check
marg_norm$solve_check
xv <- marg_norm$table[,1]
plot(xv, marg_norm$table[,2], type = "l")
plot(marg_norm$table[,1], marg_norm$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -3, to = 3, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHnorm, pars = ARCHpars)
margdens <- outer(fx_table(xvals, marg_norm$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_Gauss.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
abline(0,1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty = 3)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_norm$table)
fFxI <- FxI_table(uvals, marg_norm$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHnorm, pars = ARCHpars)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
sum(abs(copdens -1))/ngrid^2
# numerical integration is imprecise in corners
#pdf("density_ARCH_Gauss.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
abline(0,1,lty=2)
abline(v=0.5, lty=3)
abline(h=0.5, lty=3)
#dev.off()

#####################################################
# Absolute ARCH(1) normal innovations 
#####################################################

marg_absnorm <- marg_approx(c(0,5), innov = "norm", pars = ARCHpars, absolute = TRUE)
marg_absnorm$ev_check
marg_absnorm$solve_check
xv <- marg_absnorm$table[,1]
plot(xv, marg_absnorm$table[,2], type = "l")
plot(marg_absnorm$table[,1], marg_absnorm$columnsums, ylim = c(0,1))
abline(v= quantile(abs(X),c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = 0, to = 3, length=ngrid)
# Double cond prob
conddens <- 2*outer(xvals,xvals,conddensARCHnorm, pars = ARCHpars)
margdens <- outer(fx_table(xvals, marg_absnorm$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_Gauss_abs.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
abline(0,1, lty=2)
#dev.off()

# copula density picture
ngrid = 99
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_absnorm$table)
fFxI <- FxI_table(uvals, marg_absnorm$table, fx = TRUE)
# Double again
numerator <- 2*outer(FxI, FxI, conddensARCHnorm, pars = ARCHpars)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
# numerical integration is imprecise in corners
#pdf("density_ARCH_Gauss_abs.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 100)
abline(0,1,lty=2)
#dev.off()

#####################################################
# ARCH(1) with t4 innovations
#####################################################

set.seed(123)
X <- simgarch(dist = "t", n=10000, pars = ARCHpars, distpars = nuA, copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_t4 <- marg_approx(c(-5,5), innov = "t", pars = ARCHpars, innovpars = nuA)
marg_t4$ev_check
marg_t4$solve_check
xv <- marg_t4$table[,1]
plot(xv, marg_t4$table[,2], type = "l")
plot(marg_t4$table[,1], marg_t4$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -3, to = 3, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHt, pars = ARCHpars, innovpars = nuA)
margdens <- outer(fx_table(xvals, marg_t4$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_t4.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
abline(0,1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty = 3)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_t4$table)
fFxI <- FxI_table(uvals, marg_t4$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHt, pars = ARCHpars, innovpars = nuA)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
sum(abs(copdens -1))/ngrid^2
# numerical integration is imprecise in corners
#pdf("density_ARCH_t4.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
abline(0,1,lty=2)
abline(v=0.5, lty=3)
abline(h=0.5, lty=3)
#dev.off()

#####################################################
# ARCH(1) with t2.5 innovations
#####################################################

set.seed(123)
X <- simgarch(dist = "t", n=10000, pars = ARCHpars, distpars = nuB, copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_t2.5 <- marg_approx(c(-3.5,3.5), innov = "t", pars = ARCHpars, innovpars = nuB)
marg_t2.5$ev_check
marg_t2.5$solve_check
xv <- marg_t2.5$table[,1]
plot(xv, marg_t2.5$table[,2], type = "l")
plot(marg_t2.5$table[,1], marg_t2.5$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -2, to = 2, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHt, pars = ARCHpars, innovpars = nuB)
margdens <- outer(fx_table(xvals, marg_t2.5$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_t2point5.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=60)
abline(0,1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty = 3)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_t2.5$table)
fFxI <- FxI_table(uvals, marg_t2.5$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHt, pars = ARCHpars, innovpars = nuB)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
sum(abs(copdens -1))/ngrid^2
# numerical integration is imprecise in corners
#pdf("density_ARCH_t2point5.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
abline(0,1,lty=2)
abline(v=0.5, lty=3)
abline(h=0.5, lty=3)
#dev.off()

#####################################################
# ARCH(1) Gaussian innovations and leverage
#####################################################

set.seed(123)
X <- simgarch(dist = "norm", n=10000, pars = ARCHlevpars, levpar = levpar, copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_norm_lev <- marg_approx(c(-5,5), innov = "norm", pars = ARCHlevpars, levpar = levpar)
marg_norm_lev$ev_check
marg_norm_lev$solve_check
xv <- marg_norm_lev$table[,1]
plot(xv, marg_norm_lev$table[,2], type = "l")
plot(marg_norm_lev$table[,1], marg_norm_lev$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -3, to = 3, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHnorm, pars = ARCHlevpars, levpar = levpar)
margdens <- outer(fx_table(xvals, marg_norm_lev$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_Gauss_leverage.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
abline(0,1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty = 3)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_norm_lev$table)
fFxI <- FxI_table(uvals, marg_norm_lev$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHnorm, pars = ARCHlevpars, levpar =  levpar)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
# numerical integration is imprecise in corners
#pdf("density_ARCH_Gauss_leverage.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
abline(h=0.5, lty=3)
#dev.off()

#####################################################
# ARCH(1) scaled skewed Student t innovations 
#####################################################

set.seed(123)
X <- simgarch(dist = "skewt", n=10000, pars = ARCHpars, distpars = c(nuA, lambda), copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_skewt <- marg_approx(c(-5,5), innov = "skewt", pars = ARCHpars, innovpars = c(nuA, lambda))
marg_skewt$ev_check
marg_skewt$solve_check
xv <- marg_skewt$table[,1]
plot(xv, marg_skewt$table[,2], type = "l")
plot(marg_skewt$table[,1], marg_skewt$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -3, to = 3, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHskewt, pars = ARCHpars, levpar =0, innovpars =c(nuA, lambda))
margdens <- outer(fx_table(xvals, marg_skewt$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_ARCH_skewt4.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_skewt$table)
fFxI <- FxI_table(uvals, marg_skewt$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHskewt, pars = ARCHpars, innovpars = c(nuA, lambda))
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
# numerical integration is imprecise in corners
#pdf("density_ARCH_skewt4.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
#dev.off()

#####################################################
# AR(1) - ARCH(1) normal innovations 
#####################################################

set.seed(123)
X <- simgarch(dist = "norm", n=100000, pars = ARCHpars, meanpar = phi, copula = FALSE)
quantile(X, c(0.005, 0.01, 0.99, 0.995))

marg_norm_ar1 <- marg_approx(c(-5,5), innov = "norm", pars = ARCHpars, levpar = 0, meanpar = phi)
marg_norm_ar1$ev_check
marg_norm_ar1$solve_check
xv <- marg_norm_ar1$table[,1]
plot(xv, marg_norm_ar1$table[,2], type = "l")
plot(marg_norm_ar1$table[,1], marg_norm_ar1$columnsums, ylim = c(0,1))
abline(v= quantile(X,c(0.01,0.99)))

# Joint density picture
ngrid <- 99
xvals <- seq(from = -3, to = 3, length=ngrid)
conddens <- outer(xvals,xvals,conddensARCHnorm, pars = ARCHpars, levpar = 0, meanpar = phi)
margdens <- outer(fx_table(xvals, marg_norm_ar1$table), rep(1,ngrid))
jointdens <- conddens*margdens
# Joint density of (X_{t-1}, X_{t})
#pdf("xdensity_AR_ARCH_Gauss.pdf")
par(pty="s")
contour(xvals,xvals,jointdens, nlevels=30)
abline(0,1, lty=2)
#dev.off()

# copula density picture
uvals <- seq(from = 1/(ngrid+1), to = ngrid/(ngrid+1), length = ngrid)
FxI <- FxI_table(uvals, marg_norm_ar1$table)
fFxI <- FxI_table(uvals, marg_norm_ar1$table, fx = TRUE)
numerator <- outer(FxI, FxI, conddensARCHnorm, pars = ARCHpars, levpar =  0, meanpar = phi)
denominator <- outer(rep(1,ngrid), fFxI)
copdens <- numerator/denominator
# Check it approximately fulfills properties of copula density
sum(copdens)/ngrid^2
apply(copdens,1,sum)/ngrid
# numerical integration is imprecise in corners
#pdf("density_AR_ARCH_Gauss.pdf")
par(pty="s")
contour(uvals, uvals, copdens, nlevels = 50)
abline(0,1, lty=2)
#dev.off()


