library(tscopula)

# Code to generate Figure 1

mod3 <- dvinecopula(
  family = c("Clayton","Clayton", "Clayton"),
  pars = list(2, 2, 4),
  rotation = c(180, 0, 180)
)

set.seed(1)
data <- sim(mod3, n = 10000)
ts.plot(data, ylab = expression(U[t]))

# ARMA Copula Filter Example (Figure 2)

parlist <- list(ar = 0.95, ma = -0.85)

modFrank <- dvinecopula2(family = "frank",
                         pars = parlist,
                         maxlag = Inf)

modGauss <- dvinecopula2(family = "gauss",
                         pars = parlist,
                         maxlag = Inf)

modClayton <- dvinecopula2(family = "clayton",
                           pars = parlist,
                           maxlag = Inf,
                           negtau = "gauss")

modJoe <- dvinecopula2(family = "joe",
                       pars = parlist,
                       maxlag = Inf,
                       negtau = "gauss")

modGumbel <- dvinecopula2(family = "gumbel",
                          pars = parlist,
                          maxlag = Inf,
                          negtau = "gauss")

n <- 200
set.seed(17)
z <- runif(n+1)

copfilter <- function(z, mod){
  n <- length(z) -1
  U <- rep(NA, n)
  for (k in 1:n){
    pc_list <- tscopula:::mklist_dvine2(mod, k, truncate = TRUE)
    pcs <- lapply(seq_along(pc_list), function(i) {
      replicate(k - i + 1, pc_list[[i]], simplify = FALSE)
    })
    vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
    innov <- z[(length(z)-k):(length(z))]
    U[k] <- rvinecopulib::inverse_rosenblatt(t(innov), vc_short)[1,(k+1)]
  }
  U
}

UGauss <- copfilter(z, modGauss)
UFrank <- copfilter(z, modFrank)
UClayton <- copfilter(z, modClayton)
UGumbel <- copfilter(z, modGumbel)
UJoe <- copfilter(z, modJoe)

ts.plot(UGauss, ylim = range(UGauss, UFrank, UClayton, UJoe, UGumbel),
        ylab="U", xlab="k")
abline(h = UGauss[n])
lines(1:n, UFrank, col=2)
abline(h = UFrank[n], col = 2)
lines(1:n, UClayton, col=3)
abline(h = UClayton[n], col = 3)
lines(1:n, UGumbel, col=4)
abline(h = UGumbel[n], col = 4)
lines(1:n, UJoe, col=5)
abline(h = UJoe[n], col = 5)
legend(x=c(125, 200), y = c(0.7,0.83), 
       legend = c("Gauss", "Frank","Clayton","Gumbel","Joe"), col=1:5, lty = rep(1,5), bty="n")


# ARFIMA Copula Filter Example (Figure 3)
# WARNING: VERY SLOW!!!!
# Only execute if you really want to independently generate Figure 3

parlist <- list(ar = 0.95, ma = -0.85, d = 0.02)

modGauss <- dvinecopula2(family = "gauss",
                         kpacf = "kpacf_arfima",
                         pars = parlist,
                         maxlag = Inf)

modFrank <- dvinecopula2(family = "frank",
                         kpacf = "kpacf_arfima",
                         pars = parlist,
                         maxlag = Inf)

modClayton <- dvinecopula2(family = "clayton",
                           kpacf = "kpacf_arfima",
                           pars = parlist,
                           maxlag = Inf,
                           negtau = "gauss")

modJoe <- dvinecopula2(family = "joe",
                       kpacf = "kpacf_arfima",
                       pars = parlist,
                       maxlag = Inf,
                       negtau = "gauss")

modGumbel <- dvinecopula2(family = "gumbel",
                          kpacf = "kpacf_arfima",
                          pars = parlist,
                          maxlag = Inf,
                          negtau = "gauss")

n <- 700
set.seed(17)
z <- runif(n+1)
UGauss <- copfilter(z, modGauss)
UFrank <- copfilter(z, modFrank)
UClayton <- copfilter(z, modClayton)
UGumbel <- copfilter(z, modGumbel)
UJoe <- copfilter(z, modJoe)

ts.plot(UGauss, ylim = range(UGauss, UFrank, UClayton, UJoe, UGumbel),
        ylab="U", xlab="k")
abline(h = UGauss[n])
lines(1:n, UFrank, col=2)
abline(h = UFrank[n], col = 2)
lines(1:n, UClayton, col=3)
abline(h = UClayton[n], col = 3)
lines(1:n, UGumbel, col=4)
abline(h = UGumbel[n], col = 4)
lines(1:n, UJoe, col=5)
abline(h = UJoe[n], col = 5)
legend(x=c(500, 700), y = c(0.7,0.95), 
       legend = c("Gauss", "Frank","Clayton","Gumbel","Joe"), col=1:5, lty = rep(1,5), bty="n")
