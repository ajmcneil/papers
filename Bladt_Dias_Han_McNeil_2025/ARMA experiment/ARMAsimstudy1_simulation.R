library(tscopula)
library(parallel)
library(pbmcapply)
numberOfCores <- detectCores()

# selection of nice ARMA model
select_armapar <- function(tol1= 0.15, 
                           tol2 = 0.85, 
                           tol3 = 0.1, 
                           tol4 = 1e-04,
                           el = 20){
  repeat{
    phi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    psi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    tau <- kendall(armacopula(list(ar=phi, ma = psi)), lagmax = 40)
    cond1 <- (abs(phi+psi) > tol1)
    cond2 <- (min(abs(tau[1:3]))  > tol3)   
    cond3 <- (max(abs(tau[(el+1):40])) < tol4)
    if (cond1 & cond2 & cond3) break
  }
  list(ar = phi, ma = psi)
}

nsim <- 1000
mforecast <- 100
ntotal <- 1100
copchoices <- c("gumbel", "clayton", "frank", "joe", "gauss")
nreplace <- 3


## ARMA SIMULATION
set.seed(1)
system.time({
  ARMA_parameters <- mclapply(1:nsim, function(i){
    select_armapar(tol1 = 0.15, tol2 = 0.85, tol3 = 0.1, tol4 = 1e-04, el = 20)
  }, mc.cores = numberOfCores)
})

system.time({
  ARMA_Gdata <- pbmclapply(ARMA_parameters, function(pars){
    Gmod <- armacopula(pars = pars)
    qnorm(sim(Gmod, ntotal))
  }, mc.cores = numberOfCores)
})

system.time({
  ARMA_NGmods <- pbmclapply(ARMA_parameters, function(pars){
    repeat{
      copulas <- sample(copchoices, nreplace, replace = TRUE)
      if (length(copulas[copulas != "gauss"]) > 0) # make sure some are non-Gaussian
        break
    }
    posrot <- sample(c(0, 180), nreplace, replace = TRUE)
    negrot <- sample(c(90, 270), nreplace, replace = TRUE)
    sdvinecopula(family = copulas,
                          posrot = posrot,
                          negrot = negrot,
                          kpacf = "kpacf_arma",
                          pars = pars)
  }, mc.cores = numberOfCores)
})

system.time({
  ARMA_NGdata <- pbmclapply(ARMA_NGmods, function(mod){
    qnorm(sim(mod, ntotal))
  }, mc.cores = numberOfCores)
})


save(ARMA_NGmods, ARMA_NGdata, ARMA_Gdata, file = "ARMA_simstudydata.RData")


