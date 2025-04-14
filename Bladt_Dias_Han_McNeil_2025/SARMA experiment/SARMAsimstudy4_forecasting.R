library(tscopula)
library(parallel)
numberOfCores <- detectCores()
library(pbmcapply)

load("SARMA_simstudydata.RData")
load("SARMA_simstudyests.RData")

# compute quantile scores
quantile_score <- function(alpha, quants, forecast){
  weights <- ifelse(forecast <= quants, 1-alpha, alpha)
  weights*abs(forecast - quants)
}

(ndata <- as.numeric(names(SARMA_Gmods_Gest[[1]])))
J <- 20
alpha <- (1:(J-1))/J
(mforecast <- length(SARMA_NGdata[[1]]) - max(ndata))




# scores for Gaussian forecaster
system.time({
SARMA_G_scores <- pbmcmapply(function(models, alldata){
  scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                  dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
  for (j in 1:length(ndata)){
    xdata <- as.numeric(alldata[((ndata[3] - ndata[j] + 1) : ndata[3])])
    newdata <- as.numeric(alldata[(ndata[3] + 1) : (ndata[3] + mforecast)])
    for (k in 1:length(newdata)){
      udata <- strank(xdata)
      uquant <- predict(models[[j]]@tscopula, data = udata, x = alpha, type = "qf")
      xquant <- quantile(xdata, uquant)
      scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
      xdata <- c(xdata, newdata[k])[-1]
    }
  }
  scores
}, SARMA_NGmods_Gest, SARMA_NGdata, SIMPLIFY = FALSE)
})

# pits for Gaussian forecaster
system.time({
  SARMA_G_pits <- pbmcmapply(function(models, alldata){
    pits <- array(NA, dim = c(mforecast, length(ndata)),
                    dimnames = list(forecast = NULL, ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- as.numeric(alldata[((ndata[3] - ndata[j] + 1) : ndata[3])])
      newdata <- as.numeric(alldata[(ndata[3] + 1) : (ndata[3] + mforecast)])
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        pits[k,j] <- predict(models[[j]]@tscopula, data = udata, x = pedf(newdata[k], data = xdata), type = "df")
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    pits
  }, SARMA_NGmods_Gest, SARMA_NGdata, SIMPLIFY = FALSE)
})

# scores for non-Gaussian forecaster
system.time({
  SARMA_NG_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- as.numeric(alldata[((ndata[3] - ndata[j] + 1) : ndata[3])])
      newdata <- as.numeric(alldata[(ndata[3] + 1) : (ndata[3] + mforecast)])
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        uquant <- predict(models[[j]]@tscopula, data = udata, x = alpha, type = "qf")
        xquant <- quantile(xdata, uquant)
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, SARMA_NGmods_NGest, SARMA_NGdata, SIMPLIFY = FALSE)
})

# pits for non-Gaussian forecaster
system.time({
  SARMA_NG_pits <- pbmcmapply(function(models, alldata){
    pits <- array(NA, dim = c(mforecast, length(ndata)),
                  dimnames = list(forecast = NULL, ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- as.numeric(alldata[((ndata[3] - ndata[j] + 1) : ndata[3])])
      newdata <- as.numeric(alldata[(ndata[3] + 1) : (ndata[3] + mforecast)])
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        pits[k,j] <- predict(models[[j]]@tscopula, data = udata, x = pedf(newdata[k], data = xdata), type = "df")
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    pits
  }, SARMA_NGmods_NGest, SARMA_NGdata, SIMPLIFY = FALSE)
})

# Forecast score analysis

# Diebold-Mariano test
Diebold_Mariano <- function(scores1, scores2, dist = "unif"){
  if (!(identical(dim(scores1), dim(scores2))))
    stop ("Score matrices of unequal size")
  m <- dim(scores1)[1]
  alpha <- as.numeric(dimnames(scores1)[[2]])
  wts <- switch(dist,
                "none" = rep(1,length(alpha)),
                "unif" = 1/(alpha*(1-alpha)/2),
                "norm" = 1/quantile_weight(alpha))
  wtm <- matrix(wts, nrow = m, ncol = length(alpha), byrow=TRUE) 
  AWQS1i <- apply(scores1*wtm,1, mean)
  AWQS2i <- apply(scores2*wtm,1, mean)
  sigma <- sqrt(mean((AWQS1i - AWQS2i)^2))
  AWQS1 <- mean(AWQS1i)
  AWQS2 <- mean(AWQS2i)
  Tstat <- sqrt(m)*(AWQS1 - AWQS2)/sigma
  pvalue <- 2*(1-pnorm(abs(Tstat)))
  bestmethod <- paste("Method",ifelse(Tstat <= 0, 1, 2))
  list(AWQS1 = AWQS1, AWQS2 = AWQS2, sigma = sigma, Tstat=Tstat, Better = bestmethod, Null = "Forecast methods equally accurate",
       pvalue=pvalue)
}

# Compute inverse template function
quantile_weight <- function(theta, dist = "norm", inv = FALSE){
  n <- length(theta)
  output <- rep(NA, n)
  integrand <- function(x, alpha, dist){
    qfunc <- eval(parse(text = paste("q", dist, sep ="")))
    (as.numeric(qfunc(x) <= qfunc(alpha))-alpha)*(qfunc(alpha) - qfunc(x))
  }
  for (i in 1:n)
    output[i] <- integrate(integrand, lower = 0, upper = 1, 
                           alpha = theta[i], dist = dist)$value
  if (inv)
    output <- 1/output
  output
}


mforecast <- c(40,100)
AQS_results <- array(NA, dim = c(nsim, length(ndata), length(mforecast), 2, 2), 
                     dimnames = list(nsim = NULL, n = paste(ndata), m = paste(mforecast), weight = c("none","norm"), order = c("correct", "signif")))
alpha <- as.numeric(dimnames(SARMA_NG_scores)[[2]])
for (i in 1:nsim){
  for (j in (1:length(ndata))){
    NGmat <- SARMA_NG_scores[[i]][,,j]
    Gmat <- SARMA_G_scores[[i]][,,j]
    for (k in (1:length(mforecast))){
      mk<- mforecast[k]
      res_none <- Diebold_Mariano(NGmat[1:mk,], Gmat[1:mk,], dist = "none")
      res_wt <- Diebold_Mariano(NGmat[1:mk,], Gmat[1:mk,], dist = "norm")
      AQS_results[i,j,k,1,1] <- (res_none$AWQS1 < res_none$AWQS2)
      AQS_results[i,j,k,1,2] <- (res_none$AWQS1 < res_none$AWQS2) & (res_none$pvalue <= 0.05)
      AQS_results[i,j,k,2,1] <- (res_wt$AWQS1 < res_wt$AWQS2)
      AQS_results[i,j,k,2,2] <- (res_wt$AWQS1 < res_wt$AWQS2) & (res_wt$pvalue <= 0.05)
    }
  }
}

ftable(apply(100*AQS_results,c(2,3,4,5),mean), row.vars = c(2,1), col.vars = c(3,4))

# PIT analysis

mforecast <- c(40,100)
PIT_results <- array(NA, dim = c(nsim, length(ndata), length(mforecast), 2, 2), 
                     dimnames = list(nsim = NULL, n = paste(ndata), m = paste(mforecast), method = c("G","NG"), test = c("KS", "LB")))
for (i in 1:nsim){
  for (j in (1:length(ndata))){
    NGpits <- SARMA_NG_pits[[i]][,j]
    Gpits <- SARMA_G_pits[[i]][,j]
    for (k in (1:length(mforecast))){
      mk<- mforecast[k]
      PIT_results[i, j, k, "NG", "KS"] <- (ks.test(NGpits[1:mk], "punif")$p.value <= 0.05)
      PIT_results[i, j, k, "G", "KS"] <- (ks.test(Gpits[1:mk], "punif")$p.value <= 0.05)
      PIT_results[i, j, k, "NG", "LB"] <- (Box.test(NGpits[1:mk], lag = 5, type = "L")$p.value <= 0.05)
      PIT_results[i, j, k, "G", "LB"] <- (Box.test(Gpits[1:mk], lag = 5, type = "L")$p.value <= 0.05)
    }
  }
}

ftable(100*apply(PIT_results,c(2,3,4,5),mean), row.vars = c(2,1), col.vars = c(4,3))



save(SARMA_NG_scores, SARMA_NG_pits, SARMA_G_scores, SARMA_G_pits, file = "SARMA_simstudypred.RData")
  
  
