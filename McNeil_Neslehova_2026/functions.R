FZgivV <- function(z1, z2, v1, v2, copV, wtmodel){
  if (wtmodel$type == "dvine"){
    arg1 <- hbicop(cbind(v1, v2), cond_var = 1, copV)
    arg2 <- hbicop(cbind(v1, v2), cond_var = 2, copV)
    barg1 <- hbicop(cbind(z1, arg1), cond_var = 2, wtmodel$copZ1V2)
    barg2 <- hbicop(cbind(z2, arg2), cond_var = 2, wtmodel$copV1Z2)
    return(pbicop(cbind(barg1, barg2), wtmodel$copZV))
  }
  else if (wtmodel$type == "condcopula"){
    return(ifelse(wtmodel$switch(v1, v2, wtmodel$k),
                  pbicop(cbind(z1, z2), wtmodel$cop1),
                  pbicop(cbind(z1, z2), wtmodel$cop2)))
  }
  else
    stop("Unknown wtmodel type")
}

# udpf <- function(u, udp){
#   if (is(udp, "vtransform"))
#     return(udptrans(udp, u))
#   else if (is(udp, "udpcosine"))
#     return(udpfunc(u, udp$degree, udp$type))
#   else
#     stop("Type of udp function not recognized")
# }
#
# udpd <- function(u, udp){
#   if (is(udp, "Vtransform"))
#     return(tscopula::vgradient(udp, u))
#   else if ("degree" %in% names(udp))
#     return(udpderiv(u, udp$degree, udp$type))
#   else
#     stop("Type of udp function not recognized")
# }
#
# udpp <- function(udp){
#   if (is(udp, "vtransform"))
#     return(c(0, as.numeric(udp@pars)[1], 1))
#   else if ("degree" %in% names(udp)){
#     type <- udp$type
#     degree <- udp$degree
#     if (type == "cosine")
#       part <- c(0,(1:degree)/ degree)
#     else if (type == "legendre")
#       part <- c(0, Re(polyroot(basiscor:::legendre_initial(degree)$cfsD)), 1)
#     else
#       stop("Basis type not implemented")
#     return(part)
#   }
#   else
#     stop("Type of udp function not recognized")
# }
#
# bvalues <- function(v, udp, part, roottol = 1e-10, endtol = 1e-10){
#   n <- length(v)
#   if (is(udp, "VtransformI"))
#     output <- cbind(rep(0, n), as.numeric(udp@pars[1]), rep(1, n))
#   else if (is(udp, "Vtransform"))
#     output <- cbind(rep(0, n), tscopula::vdownprob(udp, v), rep(1,n))
#   else if ("degree" %in% names(udp)){
#     output <- matrix(0, nrow = n, ncol = length(part))
#     degree <- udp$degree
#     type <- udp$type
#     if (type == "legendre") {
#       y <- basiscor::qLegendre(v, degree, roottol, endtol)
#       legdata <- basiscor:::legendre_initial(degree)
#     }
#     if (type == "cosine")
#       y <- sqrt(2) * cos(pi * (1 - v))
#     for (i in 1:n){
#       rroots <- switch(type, cosine = basiscor:::cosine_roots(y[i], degree),
#                        legendre = basiscor:::legendre_roots(legdata$cfs, y[i], roottol))
#       probs <- 1/abs(udpderiv(rroots, degree, type, roottol))
#       probs <- probs/sum(probs)
#       for (j in 1:length(part)){
#         output[i,j] <- sum(probs[rroots <= part[j]])
#       }
#     }
#   }
#   else
#     stop("Type of udp function not recognized")
#   output
# }


bsicopula <- function(u1, u2, udp1, udp2, copV, wtmodel = NULL, wtonly = FALSE){
  v1 <- udptrans(udp1, u1)
  v2 <- udptrans(udp2, u2)
  maincop <- dbicop(cbind(v1, v2), family = copV)
  if (is.null(wtmodel))
    return(pmax(maincop, 0))
  u1d <- abs(udpderiv(udp1, u1))
  u2d <- abs(udpderiv(udp2, u2))
  part1 <- udp:::udpbreaks(udp1)
  part2 <- udp:::udpbreaks(udp2)
  bvalues1 <- udp:::udpbreakcdf(udp1, v1, part1)
  bvalues2 <- udp:::udpbreakcdf(udp2, v2, part2)
  Acell1 <- cut(u1, part1, labels = FALSE)
  Acell2 <- cut(u2, part2, labels = FALSE)
  b11 <- bvalues1[cbind(seq_len(nrow(bvalues1)), Acell1)]
  b12 <- bvalues1[cbind(seq_len(nrow(bvalues1)), (Acell1+1))]
  b21 <- bvalues2[cbind(seq_len(nrow(bvalues2)), Acell2)]
  b22 <- bvalues2[cbind(seq_len(nrow(bvalues2)), (Acell2+1))]
  blockp <- FZgivV(b12, b22, v1, v2, copV, wtmodel) - FZgivV(b11, b22, v1, v2, copV, wtmodel) -
    FZgivV(b12, b21, v1, v2, copV, wtmodel) + FZgivV(b11, b21, v1, v2, copV, wtmodel)
  wtfunc <- blockp * u1d * u2d
  if (wtonly)
    return(pmax(wtfunc, 0))
  pmax(maincop * wtfunc, 0)
}

negloglik_bsicopula <- function(theta, data){
  delta1 <- theta[1]
  delta2 <- theta[2]
  if ((delta1 <=0) | (delta1 >= 1) | (delta2 <= 0) | (delta2 >= 1))
    return(NA)
  udp1 <- vlinear(delta1)
  udp2 <- vlinear(delta2)
  copV <- tryCatch(bicop_dist(family = "joe", parameters = theta[3]), error = function(e) {
    return(NA)
  })
  if (is.na(copV[[1]]))
    return(NA)
  if (length(theta) == 3)
    wtmodel <- NULL
  else if (length(theta) == 4){
    copZ1V2 <- tryCatch(bicop_dist("gauss", parameters = 0), error = function(e) {
      return(NA)
    })
    copV1Z2 <- tryCatch(bicop_dist("gauss", parameters = 0), error = function(e) {
      return(NA)
    })
    copZV <- tryCatch(bicop_dist("gauss", parameters = theta[4]), error = function(e) {
      return(NA)
    })
    if (is.na(copZV[[1]]))
      return(NA)
    wtmodel <- list(type = "dvine", copZ1V2 = copZ1V2, copV1Z2 = copV1Z2, copZV = copZV)
  }
  else{
    copZ1V2 <- tryCatch(bicop_dist("gauss", parameters = theta[5]), error = function(e) {
      return(NA)
    })
    copV1Z2 <- tryCatch(bicop_dist("gauss", parameters = theta[6]), error = function(e) {
      return(NA)
    })
    copZV <- tryCatch(bicop_dist("gauss", parameters = theta[4]), error = function(e) {
      return(NA)
    })
    if (is.na(copZ1V2[[1]]) | is.na(copV1Z2[[1]]) | is.na(copZV[[1]]))
      return(NA)
    wtmodel <- list(type = "dvine", copZ1V2 = copZ1V2, copV1Z2 = copV1Z2, copZV = copZV)
  }
  sum(-log(bsicopula(data[,1], data[,2], udp1, udp2, copV, wtmodel)))
}


fit_bsicopula <- function(data){
  theta <- c(0.5, 0.5, 1.5)
  round1 <- optim(theta, negloglik_bsicopula, data = data, hessian = TRUE, control = list(maxit =1000))
  theta <- c(0.5, 0.5, round1$par[3], 0)
  round2 <- optim(theta, negloglik_bsicopula, data = data, hessian = TRUE, control = list(maxit =2000))
  theta <- c(round2$par[1:4], 0, 0)
  round3 <- optim(theta, negloglik_bsicopula, data = data, hessian = TRUE, control = list(maxit =2000))
  list(stage1 = round1$par,
       se1 = safe_ses(round1$hessian),
       stage2 = round2$par,
       se2 = safe_ses(round2$hessian),
       stage3 = round3$par,
       se3 <- safe_ses(round3$hessian),
       convergence = c(round2$convergence, round3$convergence),
       stage1_ll = -round1$value, stage2_ll = - round2$value, stage3_ll = - round3$value,
       pvalueA = 1 - pchisq(2*(round1$value - round2$value), 1),
       pvalueB = 1 - pchisq(2*(round2$value - round3$value), 2))
}

safe_ses <- function(hess) {
  hessinverse <- tryCatch(solve(hess), error = function(e) {
    warning("Hessian can't be inverted")
    return(diag(rep(NA, length(diag(hess)))))
  })
  sqrt(abs(diag(hessinverse)))
}

