library(tscopula)
library(parallel)
numberOfCores <- detectCores()
library(pbmcapply)

load("ARMA_simstudydata.RData")
load("ARMA_simstudyests.RData")

# compute quantile scores
quantile_score <- function(alpha, quants, forecast){
  weights <- ifelse(forecast <= quants, 1-alpha, alpha)
  weights*abs(forecast - quants)
}

(ndata <- as.numeric(names(ARMA_Gmods_Gest[[1]])))
J <- 20
alpha <- (1:(J-1))/J
(mforecast <- length(ARMA_NGdata[[1]]) - max(ndata))



# scores for an oracle (only 500 data) j = 2
system.time({
  ARMA_oracle_scores <- pbmcmapply(function(model, alldata){
    scores <- array(NA, dim = c(mforecast, length(alpha)),
                    dimnames = list(forecast = NULL, quantile = paste(alpha)))
      xdata <- alldata[((ndata[3] - ndata[2] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        uquant <- predict(model, data = udata, x = alpha, type = "qf")
        xquant <- quantile(xdata, uquant)
        scores[k,] <- quantile_score(alpha, xquant, newdata[k])
        xdata <- c(xdata, newdata[k])[-1]
      }
    scores
  }, ARMA_NGmods, ARMA_NGdata, SIMPLIFY = FALSE)
})

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

# pits for Gaussian forecaster
system.time({
  ARMA_G_pits <- pbmcmapply(function(models, alldata){
    pits <- array(NA, dim = c(mforecast, length(ndata)),
                    dimnames = list(forecast = NULL, ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        pits[k,j] <- predict(models[[j]]@tscopula, data = udata, x = pedf(newdata[k], data = xdata), type = "df")
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    pits
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

# pits for non-Gaussian forecaster
system.time({
  ARMA_NG_pits <- pbmcmapply(function(models, alldata){
    pits <- array(NA, dim = c(mforecast, length(ndata)),
                  dimnames = list(forecast = NULL, ndata = paste(ndata)))
    for (j in 1:length(ndata)){
      xdata <- alldata[((ndata[3] - ndata[j] + 1) : ndata[3])]
      newdata <- alldata[(ndata[3] + 1) : (ndata[3] + mforecast)]
      for (k in 1:length(newdata)){
        udata <- strank(xdata)
        pits[k,j] <- predict(models[[j]]@tscopula, data = udata, x = pedf(newdata[k], data = xdata), type = "df")
        xdata <- c(xdata, newdata[k])[-1]
      }
    }
    pits
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


# AICresults

# Order analysis
nsim <- 1000
AICresults <- array(NA, dim = c(length(ndata), 2, nsim),
                    dimnames= list(n = ndata, type = c("NG", "G"), n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    AICresults[j,1,i] <- AIC(ARMA_NGmods_NGest[[i]][[j]])
    AICresults[j,2,i] <- AIC(ARMA_NGmods_Gest[[i]][[j]])
  }
}
apply(AICresults[,1,] < AICresults[,2,],1,sum)

diffAIC <- AICresults[,2,] - AICresults[,1,]
which(diffAIC[2,]>300)

# Illustrative plot

i <- 44 # 44th sim
j <- 2 # n = 500
m <- 100
NGmat <- ARMA_NG_scores[[i]][1:m,,j]
Gmat <- ARMA_G_scores[[i]][1:m,,j]
Diebold_Mariano(NGmat, Gmat, dist = "none")
Diebold_Mariano(NGmat, Gmat, dist = "norm")

oraclemat <- ARMA_oracle_scores[[i]][1:m,]
NGave <- apply(NGmat, 2, mean)
Gave <- apply(Gmat, 2, mean)
oracleave <- apply(oraclemat, 2, mean)
# apply weighting for normal 
NGave2 <- NGave/quantile_weight(alpha, "norm")
Gave2 <- Gave/quantile_weight(alpha, "norm")
oracleave2 <- oracleave/quantile_weight(alpha, "norm")

pdf("quantile_score_plot.pdf",width=7,height=4)
par(mfrow=c(1,2),cex=0.8)
plot(alpha, NGave, type ="l",ylim = range(Gave,NGave, oracleave), ylab = "quantile score")
lines(alpha, Gave, col = "red")
lines(alpha, oracleave, col = "green")
plot(alpha, NGave2, type ="l",ylim = range(Gave2,NGave2, oracleave2), ylab = "weighted quantile score")
lines(alpha, Gave2, col = "red")
lines(alpha, oracleave2, col = "green")
dev.off()

# Forecast score analysis

mforecast <- c(40,100)
AQS_results <- array(NA, dim = c(nsim, length(ndata), length(mforecast), 2, 2), 
                     dimnames = list(nsim = NULL, n = paste(ndata), m = paste(mforecast), weight = c("none","norm"), order = c("correct", "signif")))
alpha <- as.numeric(dimnames(ARMA_NG_scores)[[2]])
for (i in 1:nsim){
  for (j in (1:length(ndata))){
    NGmat <- ARMA_NG_scores[[i]][,,j]
    Gmat <- ARMA_G_scores[[i]][,,j]
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
    NGpits <- ARMA_NG_pits[[i]][,j]
    Gpits <- ARMA_G_pits[[i]][,j]
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



save(ARMA_NG_scores, ARMA_NG_pits, ARMA_G_scores, ARMA_G_pits, ARMA_oracle_scores, file = "ARMA_simstudypred.RData")
  
  
