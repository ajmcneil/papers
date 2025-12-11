library(copula)
library(simsalapar)
library(basiscor)

findpar <- function(srho, copula, interval = c(1,10)){
  rootfunc <- function(theta, copula, srho){
    copula@parameters <- theta
    srho - rho(copula)
  }
  uniroot(rootfunc,interval, copula = copula, srho = srho)$root
}

# Parameters chosen to give Spearman rank correlation 0.7
sptrue <- 0.25
(thetaG <- findpar(sptrue, gumbelCopula(), c(1,10)))
(thetaC <- findpar(sptrue, claytonCopula(), c(0.2, 5)))
(thetaN <- findpar(sptrue, normalCopula(), c(0, 0.95)))
2*sin(pi*sptrue/6)
pars1 <- list(thetaC = thetaC, thetaG = thetaG, thetaN = thetaN)
maxN <- 6
trueC <- basiscormatrix(claytonCopula(thetaC), maxorder = maxN, symmetric = TRUE)
trueG <- basiscormatrix(claytonCopula(thetaG), maxorder = maxN, symmetric = TRUE)
trueN <- basiscormatrix(claytonCopula(thetaN), maxorder = maxN, symmetric = TRUE)
truevals1 <- list(Clayton = trueC, Gumbel = trueG, Gauss = trueN)

sptrue <- 0.75
(thetaG <- findpar(sptrue, gumbelCopula(), c(1,10)))
(thetaC <- findpar(sptrue, claytonCopula(), c(0.2, 5)))
(thetaN <- findpar(sptrue, normalCopula(), c(0, 0.95)))
2*sin(pi*sptrue/6)
pars2 <- list(thetaC = thetaC, thetaG = thetaG, thetaN = thetaN)
maxN <- 6
trueC <- basiscormatrix(claytonCopula(thetaC), maxorder = maxN, symmetric = TRUE)
trueG <- basiscormatrix(claytonCopula(thetaG), maxorder = maxN, symmetric = TRUE)
trueN <- basiscormatrix(claytonCopula(thetaN), maxorder = maxN, symmetric = TRUE)
truevals2 <- list(Clayton = trueC, Gumbel = trueG, Gauss = trueN)


pardata <- list(pars1, pars2)
truevals <- list(truevals1, truevals2)




methods <- paste("T", 1:5, sep="")
# Define sim study
varList <-
  varlist( # constructor for an object of class 'varlist' 
    n.sim = list(type="N", expr = quote(m), value = 500), 
    n = list(type="grid", value = c(20, 50, 100, 500, 1000)),
    C = list(type="grid", expr = quote(C), value = c("Clayton","Gumbel","Gauss")),
    parset = list(type="grid", value = 1:2),
    estimator = list(type="inner", expr = quote(estimator), value = methods),
    pardata = list(value = pardata),
    maxorder = list(value = maxN),
    truevals =list(value = truevals)
    )


# define sim function
doOne <- function(n, C, parset, estimator, pardata, maxorder, truevals){
  est <- rep(NA,length(estimator))
  names(est) <- estimator
  pars <- pardata[[parset]]
  cop <- switch(C,
                "Gauss" = normalCopula(pars$thetaN),
                "t" = tCopula(pars$rhoT,df=pars$df),
                "Clayton" = claytonCopula(pars$thetaC),
                "Gumbel" = gumbelCopula(pars$thetaG))
  data <- rCopula(n,cop)
  trueval <- truevals[[parset]][[C]]
  for (i in 1:length(estimator))
    est[i] <- mean(abs(basiscormatrix(data, maxorder = maxorder, symmetric = TRUE, method=estimator[i]) - trueval))
  est
}

# Check set-up 
set.seed(1)
tmp <- doOne(n=100, C="Clayton", parset = 2, estimator=c("T1","T2","T5"), pardata = pardata, maxorder = maxN, truevals=truevals)
tmp
  


# run simulation study
res <- doLapply(varList, seed="seq", doOne=doOne,  sfile="tmp.rds", monitor=TRUE)



val <- getArray(res) 

# Compute mean error
non.sim.margins <- setdiff(names(dimnames(val)), c("n.sim"))
meanabserror <- apply(val, non.sim.margins, function(v){mean(v)}) 


# Find methods with smallest meanerror
identifybest  <- function(v){
  names(v)[order(v)[1]]
}
tab1 <- ftable(apply(meanabserror,c(2,3,4),identifybest),row.vars = c("n"), col.vars=c("parset","C")) 
toLatex(tab1, 
        fontsize="normalsize",
        caption=paste("This is  the caption"),
        vlist = varList,label="table:X1")

identifybest2  <- function(v){
  v2 <- v[names(v) %in% c("T3", "T4")]
  names(v2)[order(v2)[1]]
  }
tab2 <- ftable(apply(meanabserror,c(2,3,4),identifybest2),row.vars = c("n"), col.vars=c("parset","C")) 
toLatex(tab2, 
        fontsize="normalsize",
        caption=paste("This is  the caption"),
        vlist = varList,label="table:X1")
 
# Give full numerical results for bias, sd and mae
rv <- c("C","n") # row variables
cv <- c("parset","estimator") # column variables
fmae <- ftable(meanabserror, row.vars = rv, col.vars = cv)
round(fmae,4)

# How to make a latex table
tabL <- toLatex(round(fmae,4), 
                fontsize="normalsize",
                caption=paste("This is  the caption"),
                vlist = varList,label="table:X1")


tabL




