library(tscopula)
library(xts)
library(parallel)
library(pbmcapply)
numberOfCores <- detectCores()

# selection of nice SARMA model
select_sarmapar <- function(tol1= 0.15, 
                            tol2 = 0.85, 
                            tol3 = 0.1, 
                            tol4 = 1e-04,
                            el = 20){
  repeat{
    phi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    psi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    Phi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    Psi <- runif(1, min = tol1, max = tol2)*(-1)^rbinom(1,1,0.5)
    tau <- kendall(sarmacopula(list(ar=phi, ma = psi, sar = Phi, sma = Psi), 
                               period = 4), lagmax = 40)
    cond1 <- (abs(phi+psi) > tol1)
    cond1b <- (abs(Phi+Psi) > tol1)
    cond2 <- (min(abs(tau[1:3]))  > tol3)   
    cond3 <- (max(abs(tau[(el+1):40])) < tol4)
    if (cond1 & cond1b & cond2 & cond3) break
  }
  list(ar = phi, ma = psi, sar = Phi, sma = Psi)
}

nsim <- 1000
copchoices <- c("gumbel", "clayton", "frank", "joe", "gauss")
nreplace <- 3
ntotal <- 1100 # 1000 in sample and 100 out

# make quarterly framework
yr <- rep(2001:2275, each = 4)
qtr <- rep(1:4, 225)
yrqtr <- paste(yr, qtr, sep = "-")
yrqtr <- as.yearqtr(yrqtr)

## ARMA SIMULATION
set.seed(1)
system.time({
  SARMA_parameters <- pbmclapply(1:nsim, function(i){
    select_sarmapar(tol1 = 0.15, tol2 = 0.85, tol3 = 0.1, tol4 = 1e-04, el = 20)
  }, mc.cores = numberOfCores)
})

system.time({
  SARMA_Gdata <- pbmclapply(SARMA_parameters, function(pars){
    Gmod <- sarmacopula(pars = pars, period = 4)
    data <- qnorm(sim(Gmod, ntotal))
    xts(data, yrqtr)
  }, mc.cores = numberOfCores)
})

system.time({
  SARMA_NGmods <- pbmclapply(SARMA_parameters, function(pars){
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
                          kpacf = "kpacf_sarma4",
                          pars = pars)
  }, mc.cores = numberOfCores)
})

system.time({
  SARMA_NGdata <- pbmclapply(SARMA_NGmods, function(mod){
    data <- qnorm(sim(mod, ntotal))
    xts(data, yrqtr)
  }, mc.cores = numberOfCores)
})


save(SARMA_NGmods, SARMA_NGdata, SARMA_Gdata, file = "SARMA_simstudydata.RData")


