library(tscopula)
library(parallel)
numberOfCores <- detectCores()
library(pbmcapply)
library(forecast)
load("SARMA_simstudydata.RData")


ndata <- c(200, 500, 1000)
nreplace <- 3


# ARMA ESTIMATION

# Gaussian model to Gaussian data
system.time({
  SARMA_Gmods_Gest <- pbmclapply(SARMA_Gdata, function(alldata){
    res <- vector("list", 3)
    names(res) <- paste(ndata)
    for (j in 1:length(ndata)){
      select_data <- (ndata[3] - ndata[j] + 1) : ndata[3]
      data <- alldata[select_data]
      automod <- auto.arima(data, stationary = TRUE, start.p = 1, start.q = 1)
      initmod <- tscopula::sarmavec2list(coef(automod), c(automod$arma[1:4]))
      armamod <- sarmacopula(pars = initmod, period = 4)
      res[[j]] <- fit(armamod, strank(data), control = list(maxit = 5000))
    }
    res
  }, mc.cores = numberOfCores)
})



system.time({
    SARMA_NGmods_Gest <- pbmclapply(SARMA_NGdata, function(alldata){
      res <- vector("list", 3)
      names(res) <- paste(ndata)
      for (j in 1:length(ndata)){
        select_data <- (ndata[3] - ndata[j] + 1) : ndata[3]
        data <- alldata[select_data]
        automod <- auto.arima(data, stationary = TRUE, start.p = 1, start.q = 1)
        initmod <- tscopula::sarmavec2list(coef(automod), c(automod$arma[1:4]))
        armamod <- sarmacopula(pars = initmod, period = 4)
        res[[j]] <- fit(armamod, strank(data), control = list(maxit = 5000))
      }
      res
    }, mc.cores = numberOfCores)
  })


system.time({
    SARMA_NGmods_NGest <- pbmclapply(SARMA_NGmods_Gest, function(modelset){
      res <- vector("list", 3)
      names(res) <- paste(ndata)
      for (j in 1:length(ndata)){
        order <- modelset[[j]]@tscopula@modelspec
        nr <- nreplace
        if ((order[1] > 0) & ((order[3] + order[2] + order[4]) == 0))
          nr <- min(order[1], nreplace) # low order AR case
          res[[j]] <- auto_dvine(modelset[[j]], 
                 nreplace = nr,
                 tautol = 1e-04,
                 maxlag = 25,
                 ICtol = 0.5,
                 nstrike = 3,
                 criterion = "AIC",
                 choices = c("gumbel", "clayton", "frank", "joe"),
                 verbose = FALSE)
      }
      res
    }, mc.cores = numberOfCores)
  })

# IMPORTANT RESULTS

save(SARMA_Gmods_Gest, SARMA_NGmods_Gest, SARMA_NGmods_NGest, file = "SARMA_simstudyests.RData")



