library(tscopula)
library(forecast)
load("ARMA_simstudydata.RData")
load("ARMA_simstudyests.RData")
load("ARMA_simstudyests2.RData")
load("ARMA_simstudypred.RData") 

nsim <- 1000
ndata <- c(200, 500, 1000)
nreplace <- 3




# ANALYSIS OF RESULTS

# ANALYSIS OF RESULTS

# Order analysis
orderresults <- array(NA, dim = c(length(ndata), 2, nsim),
                      dimnames= list(n = ndata, type = c("NG", "G"), n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    par <- ARMA_NGmods_Gest[[i]][[j]]@tscopula@pars
    orderresults[j,1,i] <- paste("(", paste(c(length(par$ar), length(par$ma)), collapse=","), ")", sep="")
    orderresults[j,2,i] <- paste("(", paste(as.numeric(ARMA_Gmods_Gest[[i]][[j]]@tscopula@modelspec), collapse=","), ")", sep="")
  }
}
#
apply(orderresults[,"NG",],1 , function(v,nsim){sum(v == "(1,1)")/nsim}, nsim = nsim)
apply(orderresults[,"NG",], 1, function(v,nsim){sum(v %in% c("(1,1)", "(1,0)", "(0,1)", "(1,2)", "(2,1"))/nsim}, nsim = nsim)
apply(orderresults[,"G",], 1, function(v,nsim){sum(v == "(1,1)")/nsim}, nsim = nsim)
apply(orderresults[,"G",], 1 , function(v,nsim){sum(v %in% c("(1,1)", "(1,0)", "(0,1)", "(1,2)", "(2,1"))/nsim}, nsim = nsim)



# Estimated copula analysis SEMI-PARAMETRIC
estcops <- array("none", dim = c(nreplace, length(ndata), nsim),
                 dimnames= list(cops = 1:nreplace, n = ndata, n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    mod <- ARMA_NGmods_NGest[[i]][[j]]@tscopula
    if (is(mod, "sdvinecopula")){
      family <- mod@modelspec$family
      k <- length(family)
      posrot <- mod@modelspec$posrot
      negrot <- mod@modelspec$negrot
      tau <- kendall(mod)[1:k]
      rot <- ifelse(tau < 0, negrot, posrot)
      family <- paste(family, rot, sep="")
      estcops[1:k,j,i] <- family
    }
  }
}
#
truecops <- sapply(ARMA_NGmods, FUN = function(object){
  family <- object@modelspec$family
  posrot <- object@modelspec$posrot
  negrot <- object@modelspec$negrot
  tau <- kendall(object)[1:3]
  rot <- ifelse(tau < 0, negrot, posrot)
  paste(family, rot, sep="")})
copulahits <- matrix(NA, nrow = nreplace, ncol = length(ndata))
dimnames(copulahits) <- list(paste(1:nreplace, "correct"), ndata)
for (i in 1:nreplace){
  for (j in 1:length(ndata)){
    if (i == 1)
      copulahits[i,j] <- sum(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim])/nsim
    else
      copulahits[i,j] <- sum(apply(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim], 2, all))/nsim
  }
}
copulahits




# Estimated copula analysis PARAMETRIC correct
estcops <- array("none", dim = c(nreplace, length(ndata), nsim),
                 dimnames= list(cops = 1:nreplace, n = ndata, n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    mod <- ARMA_P_est[[i]][[j]]@tscopula
    if (is(mod, "sdvinecopula")){
    family <- mod@modelspec$family
    k <- length(family)
    posrot <- mod@modelspec$posrot
    negrot <- mod@modelspec$negrot
    tau <- kendall(mod)[1:k]
    rot <- ifelse(tau < 0, negrot, posrot)
    family <- paste(family, rot, sep="")
    estcops[1:k,j,i] <- family
    }
  }
}
#
truecops <- sapply(ARMA_NGmods, FUN = function(object){
  family <- object@modelspec$family
  posrot <- object@modelspec$posrot
  negrot <- object@modelspec$negrot
  tau <- kendall(object)[1:3]
  rot <- ifelse(tau < 0, negrot, posrot)
  paste(family, rot, sep="")})
copulahits <- matrix(NA, nrow = nreplace, ncol = length(ndata))
dimnames(copulahits) <- list(paste(1:nreplace, "correct"), ndata)
for (i in 1:nreplace){
  for (j in 1:length(ndata)){
    if (i == 1)
      copulahits[i,j] <- sum(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim])/nsim
    else
      copulahits[i,j] <- sum(apply(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim], 2, all))/nsim
  }
}
copulahits

# Estimated copula analysis PARMAETRIC incorrect
estcops <- array("none", dim = c(nreplace, length(ndata), nsim),
                 dimnames= list(cops = 1:nreplace, n = ndata, n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    mod <- ARMA_Pmis_est[[i]][[j]]@tscopula
    if (is(mod, "sdvinecopula")){
      family <- mod@modelspec$family
      k <- length(family)
      posrot <- mod@modelspec$posrot
      negrot <- mod@modelspec$negrot
      tau <- kendall(mod)[1:k]
      rot <- ifelse(tau < 0, negrot, posrot)
      family <- paste(family, rot, sep="")
      estcops[1:k,j,i] <- family
    }
  }
}
#
truecops <- sapply(ARMA_NGmods, FUN = function(object){
  family <- object@modelspec$family
  posrot <- object@modelspec$posrot
  negrot <- object@modelspec$negrot
  tau <- kendall(object)[1:3]
  rot <- ifelse(tau < 0, negrot, posrot)
  paste(family, rot, sep="")})
copulahits <- matrix(NA, nrow = nreplace, ncol = length(ndata))
dimnames(copulahits) <- list(paste(1:nreplace, "correct"), ndata)
for (i in 1:nreplace){
  for (j in 1:length(ndata)){
    if (i == 1)
      copulahits[i,j] <- sum(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim])/nsim
    else
      copulahits[i,j] <- sum(apply(estcops[1:i,j,1:nsim] == truecops[1:i,1:nsim], 2, all))/nsim
  }
}
copulahits





# Parameter analysis
estimates <- array(NA, dim = c(nreplace, length(ndata), 3, nsim),
                   dimnames= list(pars = paste("tau",1:nreplace,sep=""), n = ndata, 
                                  type = c("SP", "P", "Pmis"), n.sim = NULL))
errors <- estimates
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    estimates[,j,1,i] <- kendall(ARMA_NGmods_NGest[[i]][[j]])[1:nreplace]
    estimates[,j,2,i] <- kendall(ARMA_P_est[[i]][[j]])[1:nreplace]
    estimates[,j,3,i] <- kendall(ARMA_Pmis_est[[i]][[j]])[1:nreplace]
  }
}
truetau <- sapply(ARMA_NGmods, FUN = function(object, nreplace){kendall(object)[1:nreplace]}, nreplace = nreplace)
for (j in 1:length(ndata)){
  errors[,j,1,1:nsim] <- estimates[,j,1,1:nsim] - truetau[,1:nsim]
  errors[,j,2,1:nsim] <- estimates[,j,2,1:nsim] - truetau[,1:nsim]
  errors[,j,3,1:nsim] <- estimates[,j,3,1:nsim] - truetau[,1:nsim]
}
apply(errors, c(1,2,3), function(v){sqrt(mean(v^2))})
t(apply(errors, c(2,3), function(v){sqrt(mean(v^2))}))



quantile_weight <- function(theta, dist = "norm", inv = FALSE, ...){
  n <- length(theta)
  output <- rep(NA, n)
  integrand <- function(x, alpha, dist, ...){
    qfunc <- eval(parse(text = paste("q", dist, sep ="")))
    (as.numeric(qfunc(x, ...) <= qfunc(alpha, ...))-alpha)*(qfunc(alpha, ...) - qfunc(x, ...))
  }
  for (i in 1:n)
    output[i] <- integrate(integrand, lower = 0, upper = 1, 
                           alpha = theta[i], dist = dist, ...)$value
  if (inv)
    output <- 1/output
  output
}

Diebold_Mariano <- function(scores1, scores2, dist = "unif", ...){
  if (!(identical(dim(scores1), dim(scores2))))
    stop ("Score matrices of unequal size")
  m <- dim(scores1)[1]
  alpha <- as.numeric(dimnames(scores1)[[2]])
  wts <- switch(dist,
                "none" = rep(1,length(alpha)),
                "unif" = 1/(alpha*(1-alpha)/2),
                1/quantile_weight(alpha, dist, ...))
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

nsim <- 1000
mforecast <- c(40,100)
AQS_results <- array(NA, dim = c(nsim, length(ndata), length(mforecast), 3, 2), 
                     dimnames = list(nsim = NULL, n = paste(ndata), m = paste(mforecast), 
                                     comparison = c("NG vs. G", "P vs. NG", "NG vs. P (mis)"), order = c("former chosen", "signif")))
alpha <- as.numeric(dimnames(ARMA_NG_scores)[[2]])

for (i in 1:nsim){
  for (j in (1:length(ndata))){
    NGmat <- ARMA_NG_scores[[i]][,,j]
    Gmat <- ARMA_G_scores[[i]][,,j]
    Pmat <- ARMA_P_scores[[i]][,,j]
    Pmismat <- ARMA_Pmis_scores[[i]][,,j]
    for (k in (1:length(mforecast))){
      mk<- mforecast[k]
      res_none <- Diebold_Mariano(NGmat[1:mk,], Gmat[1:mk,], dist = "none")
      res2_none <- Diebold_Mariano(Pmat[1:mk,], NGmat[1:mk,], dist = "none")
      res3_none <- Diebold_Mariano(NGmat[1:mk,], Pmismat[1:mk,], dist = "none")
      AQS_results[i,j,k,1,1] <- (res_none$AWQS1 < res_none$AWQS2)
      AQS_results[i,j,k,1,2] <- (res_none$AWQS1 < res_none$AWQS2) & (res_none$pvalue <= 0.05)
      AQS_results[i,j,k,2,1] <- (res2_none$AWQS1 < res2_none$AWQS2)
      AQS_results[i,j,k,2,2] <- (res2_none$AWQS1 < res2_none$AWQS2) & (res2_none$pvalue <= 0.05)
      AQS_results[i,j,k,3,1] <- (res3_none$AWQS1 < res3_none$AWQS2)
      AQS_results[i,j,k,3,2] <- (res3_none$AWQS1 < res3_none$AWQS2) & (res3_none$pvalue <= 0.05)
    }
  }
}

ftable(apply(100*AQS_results,c(2,3,4,5),mean), row.vars = c(2,1), col.vars = c(3,4))

for (i in 1:nsim){
  for (j in (1:length(ndata))){
    NGmat <- ARMA_NG_scores[[i]][,,j]
    Gmat <- ARMA_G_scores[[i]][,,j]
    Pmat <- ARMA_P_scores[[i]][,,j]
    Pmismat <- ARMA_Pmis_scores[[i]][,,j]
    for (k in (1:length(mforecast))){
      mk<- mforecast[k]
      res_none <- Diebold_Mariano(NGmat[1:mk,], Gmat[1:mk,], dist = "t", df=6)
      res2_none <- Diebold_Mariano(Pmat[1:mk,], NGmat[1:mk,], dist = "t", df=6)
      res3_none <- Diebold_Mariano(NGmat[1:mk,], Pmismat[1:mk,], dist = "t", df=6)
      AQS_results[i,j,k,1,1] <- (res_none$AWQS1 < res_none$AWQS2)
      AQS_results[i,j,k,1,2] <- (res_none$AWQS1 < res_none$AWQS2) & (res_none$pvalue <= 0.05)
      AQS_results[i,j,k,2,1] <- (res2_none$AWQS1 < res2_none$AWQS2)
      AQS_results[i,j,k,2,2] <- (res2_none$AWQS1 < res2_none$AWQS2) & (res2_none$pvalue <= 0.05)
      AQS_results[i,j,k,3,1] <- (res3_none$AWQS1 < res3_none$AWQS2)
      AQS_results[i,j,k,3,2] <- (res3_none$AWQS1 < res3_none$AWQS2) & (res3_none$pvalue <= 0.05)
    }
  }
}



ftable(apply(100*AQS_results,c(2,3,4,5),mean), row.vars = c(2,1), col.vars = c(3,4))








