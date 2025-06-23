library(copula)
library(parallel)
library(pbmcapply)
source("C1analysis_functions.R")

nu <- 4
lambda <- 0.8
ARCHpars <- c(0.4, 0.6, 0)
set.seed(123)
ARCHdata <- simgarch("skewt", n = 50000, pars = ARCHpars, distpars = c(nu, lambda))
GARCHpars = c(0.1, 0.3, 0.6)
set.seed(123)
GARCHdata <- simgarch("skewt", n = 50000, pars = GARCHpars, distpars = c(nu, lambda))

lev <- 0.4
ARCHpars <- c(0.4, 0.3, 0)
set.seed(123)
ARCHdatalev <- simgarch("skewt", n = 50000, pars = ARCHpars, levpar = lev, distpars = c(nu, lambda))
GARCHpars = c(0.1, 0.1, 0.6)
set.seed(123)
GARCHdatalev <- simgarch("skewt", n = 50000, pars = GARCHpars, levpar = lev, distpars = c(nu, lambda))

method_names <- c("iv_ast", "iv_Joe",
           "iv_Clayton_180", 
           "ivl_ivl_ast",  "ivl_ivl_Joe", "ivl_ivl_Clayton_180",
           "ivl_ivl_K_ast", "ivl_ivl_K_Joe", "ivl_ivl_K_Clayton_180", 
           "iv2_iv2_ast",  "iv2_iv2_Joe", "iv2_iv2_Clayton_180",
           "iv2_iv2_K_ast", "iv2_iv2_K_Joe", "iv2_iv2_K_Clayton_180")

alldata <- list(ARCHdata, ARCHdatalev, GARCHdata, GARCHdatalev)

# Results 
allresults <- pbmclapply(alldata, FUN = function(Udata){
  results <- vector("list", length(method_names)) 
  names(results) <- method_names
  for (j in 1:length(results)){
    fitfunc <- eval(parse(text = paste("fit_", method_names[j], sep="")))
    startpars <- NA
    if (method_names[j] == "ivl_ivl_K_Clayton_180")
      startpars <- c(results$ivl_ivl_Clayton_180$par.ests[1], 0, 0, results$ivl_ivl_Clayton_180$par.ests[2:3])
    if (method_names[j] == "ivl_ivl_K_Joe")
      startpars <- c(results$ivl_ivl_Joe$par.ests[1], 0, 0, results$ivl_ivl_Joe$par.ests[2:3])
    if (method_names[j] == "ivl_ivl_K_ast")
      startpars <- c(results$ivl_ivl_ast$par.ests[1], 0, 0, results$ivl_ivl_ast$par.ests[2:3])
    if (method_names[j] == "iv2_iv2_K_Clayton_180")
      startpars <- c(results$iv2_iv2_Clayton_180$par.ests[1], 0, 0, results$iv2_iv2_Clayton_180$par.ests[2:5])
    if (method_names[j] == "iv2_iv2_K_Joe")
      startpars <- c(results$iv2_iv2_Joe$par.ests[1], 0, 0, results$iv2_iv2_Joe$par.ests[2:5])
    if (method_names[j] == "iv2_iv2_K_ast")
      startpars <- c(results$iv2_iv2_ast$par.ests[1], 0, 0, results$iv2_iv2_ast$par.ests[2:5])
    results[[j]] <- fitfunc(Udata, startpars = startpars)
  }
  results
}, mc.cores = 4)

convlog <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    convlog[j,i] <- allresults[[i]][[j]]$conv
dimnames(convlog)[[1]] <- method_names
convlog

npar <- as.numeric(sapply(allresults[[1]], FUN <- function(results){results$npar}))
npar

results_asymmetry <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    results_asymmetry[j,i] <- AIC(allresults[[i]][[j]])
results_asymmetry <- cbind(npar, results_asymmetry)
dimnames(results_asymmetry) <- list(method_names, c("p"," ", "leverage", " ","leverage"))
results_asymmetry

results_asymmetry2 <- cbind(results_asymmetry[,1], apply(results_asymmetry[,-1], 2, function(v){v-min(v)}))
results_asymmetry2
xtable::xtable(results_asymmetry2, digits = 0)



