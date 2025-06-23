library(copula)
library(rvinecopulib)
library(parallel)
library(pbmcapply)
numberOfCores <- detectCores()
source("C1analysis_functions.R")

nuA  <- 4
nuB <- 2.5

ARCHpars <- c(0.4,0.6, 0)
set.seed(123)
ARCHdataG <- simgarch("norm", n = 50000, pars = ARCHpars)
set.seed(123)
ARCHdatatA <- simgarch("t", n = 50000, pars = ARCHpars, distpars = nuA)
set.seed(123)
ARCHdatatB <- simgarch("t", n = 50000, pars = ARCHpars, distpars = nuB)

GARCHpars = c(0.1, 0.3, 0.6)
set.seed(123)
GARCHdataG <- simgarch("norm", n = 50000, pars = GARCHpars)
set.seed(123)
GARCHdatatA <- simgarch("t", n = 50000, pars = GARCHpars, distpars = nuA)
set.seed(123)
GARCHdatatB <- simgarch("t", n = 50000, pars = GARCHpars, distpars = nuB)
# Collect AIC results
method_names <- c("iv_ast", "mix_t", 
           "mix_Clayton", "mix_Gumbel", "mix_Joe",
           "iv_Clayton_180", "iv_Gumbel", "iv_Joe",
           "mix_K_Clayton", "mix_K_Gumbel", "mix_K_Joe",
           "iv_K_Clayton_180", "iv_K_Gumbel", "iv_K_Joe", "iv_K_ast")

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
    if (method_names[j] == "iv_K_Gumbel")
      startpars <- c(results$iv_Gumbel$par.ests, 0,0)
    if (method_names[j] == "mix_K_Clayton")
      startpars <- c(results$mix_Clayton$par.ests, 0,0)
     if (method_names[j] == "mix_K_Joe")
       startpars <- c(results$mix_Joe$par.ests, 0,0)
    if (method_names[j] == "mix_K_Gumbel")
      startpars <- c(results$mix_Gumbel$par.ests, 0,0)
  results[[j]] <- fitfunc(Udata, startpars = startpars)
  }
  results
}, mc.cores = numberOfCores)

convlog_symmetry <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    convlog_symmetry[j,i] <- allresults[[i]][[j]]$conv
convlog_symmetry

npar <- as.numeric(sapply(allresults[[1]], FUN <- function(results){results$npar}))
npar
  
results_symmetry <- matrix(NA, nrow = length(method_names), ncol = length(alldata))
for (i in 1:length(alldata))
  for (j in 1: length(method_names))
    results_symmetry[j,i] <- AIC(allresults[[i]][[j]])
results_symmetry <-cbind(npar, results_symmetry)
dimnames(results_symmetry) <- list(method_names, c("p","Gauss", "t4", "t2.5","Gauss","t4", "t2.5"))
results_symmetry
  

results_symmetry2 <- cbind(results_symmetry[,1], apply(results_symmetry[,-1], 2, function(v){v-min(v)}))
results_symmetry2
xtable::xtable(results_symmetry2, digits = 0)

results_symmetry2["iv_Clayton_180",] - results_symmetry2["iv_K_Clayton_180",]
results_symmetry2["iv_Joe",] - results_symmetry2["iv_K_Joe",]
results_symmetry2["iv_Gumbel",] - results_symmetry2["iv_K_Gumbel",]
results_symmetry2["mix_Clayton",] - results_symmetry2["mix_K_Clayton",]
results_symmetry2["mix_Joe",] - results_symmetry2["mix_K_Joe",]
results_symmetry2["mix_Gumbel",] - results_symmetry2["mix_K_Gumbel",]

allresults[[1]]$iv_K_ast$par.ests
allresults[[2]]$iv_ast$par.ests
allresults[[3]]$iv_K_ast$par.ests
allresults[[4]]$iv_K_ast$par.ests
allresults[[5]]$iv_K_ast$par.ests
allresults[[6]]$iv_K_ast$par.ests


