library(copula)
library(parallel)
library(pbmcapply)
numberOfCores <- detectCores()
source("C1analysis_functions.R")


nuA  <- 4
nuB <- 2.5

ARCHpars <- c(0.4,0.3, 0)
lev <- 0.4
set.seed(123)
ARCHdataG <- simgarch("norm", n = 50000, pars = ARCHpars, levpar = lev)
set.seed(123)
ARCHdatatA <- simgarch("t", n = 50000, pars = ARCHpars, levpar = lev, distpars = nuA)
set.seed(123)
ARCHdatatB <- simgarch("t", n = 50000, pars = ARCHpars, levpar = lev, distpars = nuB)

0.1 + 0.6 + 0.4/2
GARCHpars = c(0.1, 0.1, 0.6)
lev <- 0.4
set.seed(123)
GARCHdataG <- simgarch("norm", n = 50000, pars = GARCHpars, levpar = lev)
set.seed(123)
GARCHdatatA <- simgarch("t", n = 50000, pars = GARCHpars, levpar = lev, distpars = nuA)
set.seed(123)
GARCHdatatB <- simgarch("t", n = 50000, pars = GARCHpars, levpar = lev, distpars = nuB)

method_names <- c("iv_ast", 
                  "iv_Joe",
           "iv_Clayton_180", 
           "iva_ast",
           "iva_Joe",
           "iva_Clayton_180", 
           "iv_K_ast", "iv_K_Joe","iv_K_Clayton_180", 
           "iva_K_ast", "iva_K_Joe", "iva_K_Clayton_180")

alldata <- list(ARCHdataG, ARCHdatatA, ARCHdatatB, GARCHdataG, GARCHdatatA, GARCHdatatB)

# Results 
allresults <- pbmclapply(alldata, FUN = function(Udata){
  results <- vector("list", length(method_names)) 
  names(results) <- method_names
  for (j in 1:length(results)){
    fitfunc <- eval(parse(text = paste("fit_", method_names[j], sep="")))
    startpars <- NA
    if (method_names[j] == "iv_K_Clayton_180")
      startpars <- c(results$iv_Clayton_180$par.ests, 0,0)
    if (method_names[j] == "iv_K_Joe")
      startpars <- c(results$iv_Joe$par.ests, 0,0)
    if (method_names[j] == "iv_K_ast")
      startpars <- c(results$iv_ast$par.ests, 0,0)
    if (method_names[j] == "iva_K_Clayton_180")
      startpars <- c(results$iva_Clayton_180$par.ests[1], 0,0, results$iva_Clayton_180$par.ests[2])
    if (method_names[j] == "iva_K_Joe")
      startpars <- c(results$iva_Joe$par.ests[1], 0,0, results$iva_Joe$par.ests[2])
    if (method_names[j] == "iva_K_ast")
      startpars <- c(results$iv_ast$par.ests[1], 0,0, results$iv_ast$par.ests[2])
    results[[j]] <- fitfunc(Udata, startpars = startpars)
  }
  results
}, mc.cores = numberOfCores)

convlog <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    convlog[j,i] <- allresults[[i]][[j]]$conv
convlog

npar <- as.numeric(sapply(allresults[[1]], FUN <- function(results){results$npar}))
npar

results_leverage <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    results_leverage[j,i] <- AIC(allresults[[i]][[j]])
results_leverage <-cbind(npar, results_leverage)
dimnames(results_leverage) <- list(method_names, c("p","Gauss", "t4", "t2.5","Gauss","t4", "t2.5"))
results_leverage

results_leverage2 <- cbind(results_leverage[,1], apply(results_leverage[,-1], 2, function(v){v-min(v)}))
results_leverage2
xtable::xtable(results_leverage2, digits = 0)


results_leverage2["iv_Clayton_180",] - results_leverage2["iv_K_Clayton_180",]
results_leverage2["iv_Joe",] - results_leverage2["iv_K_Joe",]
results_leverage2["iv_ast",] - results_leverage2["iv_K_ast",]
results_leverage2["iva_Clayton_180",] - results_leverage2["iva_K_Clayton_180",]
results_leverage2["iva_Joe",] - results_leverage2["iva_K_Joe",]
results_leverage2["iva_ast",] - results_leverage2["iva_K_ast",]

