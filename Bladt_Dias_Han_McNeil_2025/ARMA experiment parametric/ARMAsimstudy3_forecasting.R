library(tscopula)
library(parallel)
numberOfCores <- detectCores()
library(pbmcapply)

load("ARMA_simstudydata.RData")
load("ARMA_simstudyests.RData")
load("ARMA_simstudyests2.RData")

# compute quantile scores
quantile_score <- function(alpha, quants, forecast){
  weights <- ifelse(forecast <= quants, 1-alpha, alpha)
  weights*abs(forecast - quants)
}

(ndata <- as.numeric(names(ARMA_Gmods_Gest[[1]])))
J <- 20
alpha <- (1:(J-1))/J
(mforecast <- length(ARMA_NGdata[[1]]) - max(ndata))

# scores for Gaussian forecaster
system.time({
  ARMA_G_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        uquant <- predict(models[[j]]@tscopula, data = udata, x = alpha, type = "qf")
        xquant <- quantile(xdata, uquant)
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, ARMA_NGmods_Gest, ARMA_NGdata, SIMPLIFY = FALSE)
})


# scores for non-Gaussian forecaster
system.time({
  ARMA_NG_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        uquant <- predict(models[[j]]@tscopula, data = udata, x = alpha, type = "qf")
        xquant <- quantile(xdata, uquant)
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, ARMA_NGmods_NGest, ARMA_NGdata, SIMPLIFY = FALSE)
})

# scores for correct parametric forecaster
system.time({
  ARMA_P_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        xquant <- predict(models[[j]], data = xdata, x = alpha, type = "qf")
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, ARMA_P_est, ARMA_NGdata, SIMPLIFY = FALSE)
})


# scores for incorrect parametric forecaster
system.time({
  ARMA_Pmis_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        xdata[pmarg(models[[j]]@margin, xdata) == 0] <- min(xdata[pmarg(models[[j]]@margin, xdata) > 0])
        xdata[pmarg(models[[j]]@margin, xdata) == 1] <- max(xdata[pmarg(models[[j]]@margin, xdata) < 1])
        xquant <- predict(models[[j]], data = xdata, x = alpha, type = "qf")
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, ARMA_Pmis_est, ARMA_NGdata, SIMPLIFY = FALSE)
})


save(ARMA_G_scores, ARMA_NG_scores, ARMA_P_scores, ARMA_Pmis_scores, file = "ARMA_simstudypred.RData")


  
