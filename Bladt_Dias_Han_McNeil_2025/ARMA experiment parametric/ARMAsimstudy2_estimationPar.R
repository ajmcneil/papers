library(tscopula)
library(parallel)
numberOfCores <- detectCores() -1
library(pbmcapply)
library(forecast)
load("ARMA_simstudydata.RData")


ndata <- c(200, 500, 1000)
nreplace <- 3


# ARMA ESTIMATION

system.time({
    ARMA_P_est <- pbmclapply(ARMA_NGdata, function(alldata){
      res <- vector("list", 3)
      names(res) <- paste(ndata)
      for (j in 1:length(ndata)){
        select_data <- (ndata[3] - ndata[j] + 1) : ndata[3]
        data <- alldata[select_data]
        Ndata <- qnorm(strank(data))
        automod <- auto.arima(Ndata, stationary = TRUE, start.p = 1, start.q = 1)
        initmod <- tscopula::armavec2list(coef(automod), automod$arma[1:2])
        armamod <- armacopula(pars = initmod)
        margmod <- fit(margin("st"), data)
        U <- pmarg(margmod, data)
        copmod <- fit(armamod, U, control = list(maxit = 5000))
        order <- copmod@tscopula@modelspec
        nr <- nreplace
        if ((order[1] > 0) & (order[2] == 0))
          nr <- min(order[1], nreplace) # low order AR case
        newcopmod <- auto_dvine(copmod, 
                                nreplace = nr,
                                tautol = 1e-04,
                                maxlag = 25,
                                ICtol = 0.5,
                                nstrike = 3,
                                criterion = "AIC",
                                choices = c("gumbel", "clayton", "frank", "joe"),
                                verbose = FALSE)
        res[[j]] <- tscm(newcopmod, margmod)
      }
      res
    }, mc.cores = numberOfCores)
  })

system.time({
  ARMA_Pmis_est <- pbmclapply(ARMA_NGdata, function(alldata){
    res <- vector("list", 3)
    names(res) <- paste(ndata)
    for (j in 1:length(ndata)){
      select_data <- (ndata[3] - ndata[j] + 1) : ndata[3]
      data <- alldata[select_data]
      Ndata <- qnorm(strank(data))
      automod <- auto.arima(Ndata, stationary = TRUE, start.p = 1, start.q = 1)
      initmod <- tscopula::armavec2list(coef(automod), automod$arma[1:2])
      armamod <- armacopula(pars = initmod)
      margmod <- fit(margin("norm"), data)
      U <- pmarg(margmod, data)
      U[U <= 1e-10] <- 1e-10
      U[U >= 1 - 1e-10] <- 1 - 1e-10
      copmod <- fit(armamod, U, control = list(maxit = 5000))
      order <- copmod@tscopula@modelspec
      nr <- nreplace
      if ((order[1] > 0) & (order[2] == 0))
        nr <- min(order[1], nreplace) # low order AR case
      newcopmod <- auto_dvine(copmod, 
                              nreplace = nr,
                              tautol = 1e-04,
                              maxlag = 25,
                              ICtol = 0.5,
                              nstrike = 3,
                              criterion = "AIC",
                              choices = c("gumbel", "clayton", "frank", "joe"),
                              verbose = FALSE)
      res[[j]] <- tscm(newcopmod, margmod)
    }
    res
  }, mc.cores = numberOfCores)
})


# IMPORTANT RESULTS

save(ARMA_P_est, ARMA_Pmis_est, file = "ARMA_simstudyests2.RData")



