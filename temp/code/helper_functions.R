


IPS <- function(rho, A, eps=1e-8){
  n <- nrow(A)
  diag(A) <- 0
  T <- diag(n) + rho*A

  W <- clique_part(A)
  Q <- IPS2(rho=rho, n=n, T=T, W=W, maxit=100, eps=eps)$Q

  return(Q)   
}


clique_part <- function(A){
  A_NELO <- as(A, "graphNEL")
  Clist  <- getCliques(A_NELO)
  Clist  <- lapply(Clist, function(x) sort(as.numeric(x)))

  V      <- sum(lower.tri(A)*A)
  NN     <- length(Clist)
  Nc     <- NN/(1:NN) 
  NcTemp <- Nc[which(Nc %% 1 == 0)]
  M      <- NcTemp[-c(1, length(NcTemp))]

  Mcomp  <- sapply(M, clique_select,Clist=Clist, V=V, NN=NN)

  nc     <- M[which(Mcomp==min(Mcomp))]
  nC     <- length(Clist)/nc
  W      <- lapply(1:nc, function(i) Clist[(i*nC-nC+1):(i*nC)])

  #Need to adjust index since C++ index begins with zero
  W1  <- lapply(1:nc, function(i) lapply(W[[i]], function(x) x-1))  
  W1
}


clique_select <- function(m, Clist, V, NN){
  nC  <- length(Clist)/m
  Wt  <- lapply(1:m, function(i) Clist[(i*nC-nC+1):(i*nC)])
  Umstar <- max(sapply(1:m, function(j) sum(sapply(Wt[[j]], length))))
  log( m*V^3 + NN * Umstar^3 )
}


negllk_nu <- function(nu, Y, X, offset, beta){
  mu   <- exp(offset + X%*%beta)

  e1 <- nu * mu * log(nu) - lgamma(nu*mu) + lgamma(Y + nu * mu) - (Y + nu*mu) * log(1+nu)
  -sum(e1)
}


negllk_eta <- function(pars, Y, X, offset){
  nu   <- exp(pars[1])
  beta <- pars[-1]
  mu   <- exp(offset + X%*%beta)

  e1 <- nu * mu * log(nu) - lgamma(nu*mu) + lgamma(Y + nu * mu) - (Y + nu*mu) * log(1+nu)
  -sum(e1)
}


negllk_rho <- function(rho, Y, X, A, W, offset, nu, beta, lambda){
  Q    <- IPS(rho=rho, A=A)
  detQ <- determinant(Q, logarithm=TRUE)$modulus[1]

  n      <- nrow(X)
  mu     <- exp(offset + X%*%beta)
  U      <- pgamma(lambda, shape=nu*mu, rate=nu)
  z      <- qnorm(U)
    z    <- ifelse(z == -Inf, -10000, z)
  qform  <- as.numeric(t(z)%*%Q%*%z)
  e1     <- 0.5 * detQ - 0.5*qform 
  -e1
}

negllk_lambda <- function(lambda, Y, X, offset, nu, beta){
  mu     <- exp(offset + X%*%beta)
  e1     <- sum(dgamma(lambda, shape=nu*mu, rate=nu, log=TRUE))
  e2     <- sum(dpois(Y,lambda,log=TRUE))
  -sum(e1 + e2)
}


negllk_hess <- function(pars, Y, X, A, W, offset){
  nu     <- exp(pars[1])
  rho    <- pnorm(pars[2])
  beta   <- pars[-c(1:2)]
  Q      <- IPS(rho=rho, A=A)
  detQ   <- determinant(Q, logarithm=TRUE)$modulus[1]

  n      <- nrow(X)
  p      <- ncol(X)
  mu     <- exp(offset + X%*%beta)
  lambda <- sapply(1:N, function(i){
                          optimize(negllk_lambda, c(10e-12, 10^4), 
                                   Y=Y[i], 
                                   X=X[i,], 
                                   offset=offset[i], 
                                   nu=nu, 
                                   beta=beta)$min})
  U      <- pgamma(lambda, shape=nu*mu, rate=nu)
  z      <- qnorm(U)
    z    <- ifelse(z == -Inf, -10000, z)
  qform  <- t(z) %*% ( Q - diag(n) ) %*% z
  e1     <- 0.5 * detQ - 0.5*qform
  e2     <- sum(dgamma(lambda, shape=nu*mu, rate=nu, log=TRUE))
  e3     <- sum(dpois(Y,lambda,log=TRUE))
  -sum(e1 + e2 + e3)
}


boot.sim.helper <- function(X, offset, beta, nu, rootQ){
  N     <- nrow(X)
  mu    <- exp(offset + X%*%beta)

  z      <- rnorm(N)
  Z      <- solve(rootQ, z)
  U      <- sapply(Z, pnorm)             #G( lambda )
  lambda <- qgamma(U, shape=mu*nu, rate=nu)
  Y      <- rpois(N, lambda) 
  Y
}


copCS <- function(Y, X, offset, A, conf.level=0.95, conf.int=0, boot.iter=1){

  W <- clique_part(A)

  ###Calculate initial values for beta and nu
  beta0  <- summary(glm(Y~X-1,family=poisson, offset=offset))$coef[,1]
  nu0    <- optimize(negllk_nu, c(0.001,10), 
                     Y=Y, X=X, offset=offset, beta=beta0)$minimum

  ###Calculate parameters
  pars    <- c(nu0,beta0)
  etafit  <- optim(pars, negllk_eta, Y=Y, X=X, offset=offset)
  nu2     <- exp(etafit$par[1])
  beta2   <- etafit$par[-1]
  mu2     <- exp(log(E) + X%*%beta2)
  lambda2 <- sapply(1:N, function(i){
                           optimize(negllk_lambda, c(10e-12, 10^4), 
                                    Y=Y[i], 
                                    X=X[i,], 
                                    offset=offset[i], 
                                    nu=nu2, 
                                    beta=beta2)$min})
  rhofit <- optimize(negllk_rho, interval=c(0.01,0.99), 
                     Y=Y, 
                     X=X, 
                     A=A, 
                     W=W, 
                     offset=offset, 
                     nu=nu2, 
                     beta=beta2, 
                     lambda=lambda2)
  rho2   <- rhofit$minimum
  pars   <- c(log(nu2), qnorm(rho2), beta2)

  ff <- optim(pars, negllk_hess, Y=Y, X=X, A=A, W=W, offset=offset)

  if(conf.int==0){
    pars[1] <- exp(pars[1])
    pars[2] <- pnorm(pars[2])
    names(pars) <- c("nu", "rho", colnames(X))
    return(pars)
  } else if(conf.int==1){
    conf.alpha <- qnorm(1-(1-conf.level)/2)
  
    ###Calculate inverse hessian (the bread)
    hess      <- hessian(negllk_hess, pars, Y=Y, X=X, A=A, W=W, offset=offset)
    I.hat.inv <- solve(hess)
   
    ###Calculate Jhat matrix (the ham)
    npar     <- length(pars)
    Q        <- IPS(rho=rho2, A=A)
    rootQ    <- chol(Q)
    J.hat    <- matrix(0,npar,npar)
    b        <- boot.iter
 
    for(j in 1:b){
      Y_    <- boot.sim.helper(X=X, offset=offset, beta=beta2, nu=nu2, rootQ=rootQ)
      gr    <- - grad(negllk_hess, pars, Y=Y_, X=X, A=A, W=W, offset=offset, method="simple")
      J.hat <- J.hat + gr %o% gr / b
    }

    G.hat.inv <- I.hat.inv %*% J.hat %*% I.hat.inv
    SDs       <- sqrt( diag( G.hat.inv ) )
    CIs       <- t(sapply(1:npar, function(i) pars[i] + c(-1,1) * conf.alpha * SDs[i]))
    res       <- cbind(pars, SDs, CIs)   

    ### Format results
    res[1,c(1,3:4)] <- exp(res[1,c(1,3:4)])
    res[2,c(1,3:4)] <- pnorm(res[2,c(1,3:4)])
    row.names(res)  <- NULL
    res             <- round(res,3)
    res             <- as.data.frame(res)
    res[,3]         <- paste0("(",res[,3],", ",res[,4],")")
    res             <- res[,-4]
    colnames(res)   <- c("Estimate", "Std.Error", 
                         paste0(100*conf.level, "% CI"))
    rownames(res)   <- c("nu", "rho", paste0("X",colnames(X)))
  }
  
 
  return(res)
}
