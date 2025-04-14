library(tscopula)
library(forecast)
load("ARMA_simstudydata.RData")
load("ARMA_simstudyests.RData")

nsim <- 1000
ndata <- c(200, 500, 1000)
nreplace <- 3




# ANALYSIS OF RESULTS

# Order analysis
orderresults <- array(NA, dim = c(length(ndata), 2, nsim),
                      dimnames= list(n = ndata, type = c("NG", "G"), n.sim = NULL))
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    par <- ARMA_NGmods_NGest[[i]][[j]]@tscopula@pars
    orderresults[j,1,i] <- paste("(", paste(c(length(par$ar), length(par$ma)), collapse=","), ")", sep="")
    orderresults[j,2,i] <- paste("(", paste(as.numeric(ARMA_Gmods_Gest[[i]][[j]]@tscopula@modelspec), collapse=","), ")", sep="")
    }
}
#
apply(orderresults[,"NG",],1 , function(v,nsim){sum(v == "(1,1)")/nsim}, nsim = nsim)
apply(orderresults[,"NG",], 1, function(v,nsim){sum(v %in% c("(1,1)", "(1,0)", "(0,1)", "(1,2)", "(2,1"))/nsim}, nsim = nsim)
apply(orderresults[,"G",], 1, function(v,nsim){sum(v == "(1,1)")/nsim}, nsim = nsim)
apply(orderresults[,"G",], 1 , function(v,nsim){sum(v %in% c("(1,1)", "(1,0)", "(0,1)", "(1,2)", "(2,1"))/nsim}, nsim = nsim)



# Estimated copula analysis
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



# Parameter analysis
estimates <- array(NA, dim = c(nreplace, length(ndata), 3, nsim),
                   dimnames= list(pars = paste("tau",1:nreplace,sep=""), n = ndata, 
                                  type = c("G-G", "NG-G", "NG-NG"), n.sim = NULL))
errors <- estimates
for (i in 1:nsim){
  for (j in 1:length(ndata)){
    estimates[,j,1,i] <- kendall(ARMA_Gmods_Gest[[i]][[j]])[1:nreplace]
    estimates[,j,2,i] <- kendall(ARMA_NGmods_Gest[[i]][[j]])[1:nreplace]
    estimates[,j,3,i] <- kendall(ARMA_NGmods_NGest[[i]][[j]])[1:nreplace]
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





# vl <- varlist(
#   n.sim = list(value = nsim), # type = N
#   pars  = list(type="grid", value = paste("tau",1:nreplace,sep="")),
#   n     = list(type="grid", value = ndata), # sample sizes
#   type = list(type="grid", value = c("NG","G")))
# 
# pdf("tmp2.pdf")
# par(cex=0.8)
# mayplot(errors,  vl, row.vars ="pars", col.vars = "type", xvar ="n")
# dev.off()







