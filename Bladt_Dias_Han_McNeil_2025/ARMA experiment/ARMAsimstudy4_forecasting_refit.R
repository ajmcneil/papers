library(tscopula)
library(parallel)
numberOfCores <- detectCores()
library(pbmcapply)

load("ARMA_simstudydata.RData")
load("ARMA_simstudyests.RData")
load("ARMA_simstudypred.RData")

# compute quantile scores
quantile_score <- function(alpha, quants, forecast){
  weights <- ifelse(forecast <= quants, 1-alpha, alpha)
  weights*abs(forecast - quants)
}

(ndata <- as.numeric(names(ARMA_Gmods_Gest[[1]])))
J <- 20
alpha <- (1:(J-1))/J
(mforecast <- length(ARMA_NGdata[[1]]) - max(ndata))

# scores for non-Gaussian forecaster with refit
system.time({
  ARMA_NG_refit_scores <- pbmcmapply(function(models, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha), length(ndata)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha), ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        model <- models[[j]]@tscopula
        if (k %% 5 == 0)
          model <- fit(model, udata)@tscopula
        uquant <- predict(model, data = udata, x = alpha, type = "qf")
        xquant <- quantile(xdata, uquant)
        scores[k,,j] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    scores
  }, ARMA_NGmods_NGest, ARMA_NGdata, SIMPLIFY = FALSE)
})



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



# Forecast score analysis

nsim <- 1000
mforecast <- c(40,100)
AQS_results <- array(NA, dim = c(nsim, length(ndata), length(mforecast), 2, 2), 
                     dimnames = list(nsim = NULL, n = paste(ndata), m = paste(mforecast), weight = c("none","norm"), order = c("refit better", "signif")))
alpha <- as.numeric(dimnames(ARMA_NG_scores)[[2]])
for (i in 1:nsim){
  for (j in (1:length(ndata))){
    norefit <- ARMA_NG_scores[[i]][,,j]
    refit <- ARMA_NG_refit_scores[[i]][,,j]
    for (k in (1:length(mforecast))){
      mk<- mforecast[k]
      res_none <- Diebold_Mariano(refit[1:mk,], norefit[1:mk,], dist = "none")
      res_wt <- Diebold_Mariano(refit[1:mk,], norefit[1:mk,], dist = "norm")
      AQS_results[i,j,k,1,1] <- (res_none$AWQS1 < res_none$AWQS2)
      AQS_results[i,j,k,1,2] <- (res_none$AWQS1 < res_none$AWQS2) & (res_none$pvalue <= 0.05)
      AQS_results[i,j,k,2,1] <- (res_wt$AWQS1 < res_wt$AWQS2)
      AQS_results[i,j,k,2,2] <- (res_wt$AWQS1 < res_wt$AWQS2) & (res_wt$pvalue <= 0.05)
    }
  }
}

ftable(apply(100*AQS_results,c(2,3,4,5),mean), row.vars = c(2,1), col.vars = c(3,4))


  
  
