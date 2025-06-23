#GARCH/ARCH simulation functions

GARCHvolfunc <- function(x, s, pars, levpar = 0){
  sqrt(pars[1] + (pars[2] + levpar*as.numeric(x < 0))*x^2 + pars[3]*s^2)
}
#
simgarch <- function(dist, n, pars, levpar = 0, distpars = NULL, 
                     meanpar = 0, copula = TRUE, outputsigma = FALSE){   
  Z <- switch(dist,
              norm = rnorm(n),
              t = {
                df <- distpars[1]
                rt(n, df)*sqrt((df-2)/df)
              },
              skewt = {
                df <- distpars[1]
                gamma <- distpars[2]
                M1 <- dt(0,df=df)*df/(df-1)*2 
                dist.mn <- M1*(gamma-1/gamma)
                M2 <- df/(df-2)
                dist.var <- (M2-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - M2
                unscaled <- skewt::rskt(n, df = df, gamma= gamma)
                (unscaled-dist.mn)/sqrt(dist.var)
              }
  )
  sigma <- rep(NA, n)
  X <- rep(NA, n)
  sigma[1] <- 1
  X[1] <- Z[1]*sigma[1]
  for (i in 2:n){
    sigma[i] <- GARCHvolfunc(X[i-1], sigma[i-1], pars, levpar) 
    X[i] <- meanpar*X[i-1] + sigma[i] * Z[i]
  }
  if (outputsigma)
    X <- sigma
  if (copula)
    output <- tscopula::strank(X)
  else 
    output <- X
  output
}

marg_approx <- function(limits, innov = "norm", pars = ARCHpars, levpar = 0, meanpar = 0, innovpars = NA,
                        granularity = 2000, absolute = FALSE){
  #  browser()
  delta <- (limits[2] - limits[1])/granularity
  xv <- seq(from = limits[1], to = limits[2], length = granularity)
  conddens <- switch(innov,
                     "norm" = conddensARCHnorm,
                     "t" = conddensARCHt,
                     "skewt" = conddensARCHskewt)
  At <- outer(xv, xv, conddens, pars = pars, innovpars = innovpars, levpar = levpar, meanpar = meanpar)
  A <- t(At)*delta
  if (absolute)
    A <- 2*A
  # column sums should be close to 1
  tmp <- eigen(A)
  evalues <- tmp$values
  # max eigenvalue should be close to 1
  ev_check <- c(max(abs(evalues)), evalues[1])
  vs <- tmp$vectors
  fx <- rep(NA,  granularity)
  for (i in 1:granularity){
    fx[i] <- -Re(vs[i,1])
  }
  fx <- fx/(sum(fx)*delta)
  solve_check <- max(A %*% fx -fx)
  list(table = cbind(x=xv,fx = fx, Fx = cumsum(fx)*delta), ev_check = ev_check, solve_check = solve_check, columnsums = apply(A,2,sum))
}

conddensARCHnorm <- function(x,y, pars, levpar = 0, meanpar = 0, innovpars = NA){
  sigV <- GARCHvolfunc(x, 0, pars, levpar = levpar)
  arg <- y - meanpar*x
  dnorm(arg/sigV)/sigV
}

conddensARCHt <- function(x,y, pars, levpar = 0, meanpar = 0, innovpars){
  nu <- innovpars[1]
  scale <- sqrt((nu-2)/nu)
  sigV <- GARCHvolfunc(x, 0, pars, levpar = levpar)
  arg <- y - meanpar*x
  dt(arg/(sigV*scale), df = nu)/(sigV*scale)
}

conddensARCHskewt <- function(x,y, pars, levpar = 0, meanpar = 0, innovpars){
  nu <- innovpars[1]
  lambda <- innovpars[2]
  M1 <- dt(0,df=nu)*nu/(nu-1)*2 
  dist.mn <- M1*(lambda-1/lambda)
  M2 <- nu/(nu-2)
  dist.sd <- sqrt((M2-M1^2)*((lambda^2)+1/(lambda^2)) + 2*(M1^2) - M2)
  sigV <- GARCHvolfunc(x, 0, pars, levpar = levpar)
  arg <- y - meanpar*x
  skewt::dskt(dist.mn + dist.sd*arg/sigV, df = nu, gamma = lambda)*dist.sd/sigV
}

fx_table <- function(x, table){
  idx <- rep(NA, length(x))
  for (j in 1:length(x)){
    idx[j] <- which.min(abs(x[j] - table[,"x"]))
  }
  table[idx,"fx"]
}

FxI_table <- function(u, table, fx = FALSE){
  idx <- rep(NA, length(u))
  for (j in 1:length(u)){
    idx[j] <- which.min(abs(u[j] - table[,"Fx"]))
  }
  if (fx)
    out <- table[idx, "fx"]
  else
    out <- table[idx, "x"]
  out
}