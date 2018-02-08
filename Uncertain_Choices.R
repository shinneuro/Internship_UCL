time_started <- sprintf('Started at: %s',Sys.time())

real_dat <- read.csv("Speekenbrink_Konstantinidis_2016_Exp1_anonymous.csv")
alldat <- subset(real_dat,id2 != 1) # excluding first participant. reindexing needed.

LinearUtil <- function(R) {
  R
}

ProspectUtil <- function(alpha,lambda,R) {
  ut <- vector("double",length=length(R))
  ut[R < 0] <- -lambda*abs(R[R<0])^alpha # only for the elements in vector R that are negative.
  ut[R >= 0] <- R[R >= 0 ]^alpha
  ut
}

bayesUpdate <- function(sg2_zeta,sg2_epsilon,id,trial,choice,util){
  n <- length(choice)
  E <- matrix(0.0,nrow=n,ncol=4)
  K <- matrix(0.0,nrow=n,ncol=4)
  S <- matrix(0.0,nrow=n,ncol=4)
  
  # E[0]=0, S[0]=1000
  E[1,] <- c(0, 0, 0, 0)
  S[1,] <- c(1000, 1000, 1000, 1000)
  
  for(i in 1:(n-1)) {
    K[i,] <- (S[i,]+sg2_zeta) / (S[i,]+sg2_zeta+sg2_epsilon)
    E[i+1,] <- E[i,] + (choice[i]==1:4)*K[i,]*(util[i] - E[i,])
    S[i+1,] <- (1-(choice[i]==1:4)*K[i,])*(S[i,]+sg2_zeta)
  }
  return(list(E,S))
}


# Probability of Maximum Utility ------------------------------------------

library(OpenMx)

PMU <- function(E,S, sg2_epsilon) {

  # fkf package used in Speekenbrink lab.
  
  # Matrix A
  
  c0 <- c(1,1,1)
  c1 <- c(-1,0,0)
  c2 <- c(0,-1,0)
  c3 <- c(0,0,-1)
  
  A <- vector("list", length = 4)
  A[[1]] <- cbind(c0,c1,c2,c3)
  A[[2]] <- cbind(c1,c0,c2,c3)
  A[[3]] <- cbind(c1,c2,c0,c3)
  A[[4]] <- cbind(c1,c2,c3,c0)
  
  p <- matrix(0,nrow=200,ncol=4)
  
  # mean and covariance
  for(i in 1:200) { # of trials
  	for( j in 1:4) { # of arms
  	  means <- as.vector(A[[j]] %*% E[i,])
  		cov <- A[[j]]%*%diag(S[i,]+sg2_epsilon)%*%t(A[[j]])
  		p[i,j] <- omxMnor(cov,means,c(0,0,0),c(Inf,Inf,Inf)) # or try pmvnorm from the mvtnorm
  		
  	}
  }
  return(p)
}

P_BayesPMU_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  alpha <- exp(par[3])
  lambda <- exp(par[4])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- PMU(E,S,sg2_epsilon)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}



# Probability of Improvement ----------------------------------------------

PI <- function(E,S) {
  p <- pnorm((E - apply(E,1,max))/sqrt(S))
  return(p)
}

P_BayesPI_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  alpha <- exp(par[3])
  lambda <- exp(par[4])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- PI(E,S)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_BayesPI_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  ut <- LinearUtil(data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- PI(E,S)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}


# Upper Confidence Bound --------------------------------------------------

UCB <- function(E,S,theta) {
  p <- exp(theta*(E+1.96*sqrt(S)))
  p <- p/rowSums(p)
  return(p)
}

P_BayesUCB_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  theta <- exp(par[3])
  alpha <- exp(par[4])
  lambda <- exp(par[5])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- UCB(E,S,theta)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_BayesUCB_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  ut <- LinearUtil(data$payoff)
  theta <- exp(par[3])
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- UCB(E,S,theta)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}


# Expected Improvement ----------------------------------------------------

EI <- function(E,S,theta) {
  m <- (E - apply(E,1,max))
  p <- exp(theta*(m*pnorm(m/sqrt(S))+sqrt(S)*dnorm(m/sqrt(S))))
  p <- p/rowSums(p)
  return(p)
}


P_BayesEI_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  theta <- exp(par[3])
  alpha <- exp(par[4])
  lambda <- exp(par[5])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- EI(E,S,theta)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_BayesEI_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  theta <- exp(par[3])
  ut <- LinearUtil(data$payoff)
  res <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  E <- res[[1]]
  S <- res[[2]]
  p <- EI(E,S,theta)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}



# Optimization ------------------------------------------------------------


P_kalman_optfun <- function(x,llfun,lower,upper) { # for each participant, parameter estimation
  nstart <- 1000
  dat <- subset(alldat,id2==x)
  
  # pars: sg_zeta, sg_epsilon, alpha, lambda
  pars <- t(runif.sobol(nstart,length(lower),1,0,123345)) # t: matrix transpose # runif: generates random deviates # runif.sobol(n, dimension, init, scrambling, seed): Uniform scrambled Sobol sequence
  pars <- t(lower + (upper - lower)*pars)
  # pars[,1:2]: sg_zeta, sg_epsilon
  pars[,1:2] <- log(pars[,1:2])
  # pars[,3:4]: alpha, lambda
  pars[,3:ncol(pars)] <- log(pars[,3:ncol(pars)])
  
  tll <- rep(NA,nstart)
  for(i in 1:nstart) {
    tll[i] <-  llfun(par=pars[i,],data=dat)
  }
  par <- pars[which.min(tll),]
  out <- optim(par=par,fn=llfun,data=dat,control=list(maxit=1000)) # optim: general-purpose optimization. Nelder-Mead simplex algorithm.
  cat("Subject",x,"estimated \n")
  return(out)
}

L_kalman_optfun <- function(x,llfun,lower,upper) { # for each participant, parameter estimation
  nstart <- 1000
  dat <- subset(alldat,id2==x)
  
  # pars: sg_zeta, sg_epsilon, alpha, lambda
  pars <- t(runif.sobol(nstart,length(lower),1,0,123345)) # t: matrix transpose # runif: generates random deviates # runif.sobol(n, dimension, init, scrambling, seed): Uniform scrambled Sobol sequence
  pars <- t(lower + (upper - lower)*pars)
  # pars[,1:2]: sg_zeta, sg_epsilon
  pars[,1:2] <- log(pars[,1:2])

  tll <- rep(NA,nstart)
  for(i in 1:nstart) {
    tll[i] <-  llfun(par=pars[i,],data=dat)
  }
  par <- pars[which.min(tll),]
  out <- optim(par=par,fn=llfun,data=dat,control=list(maxit=1000)) # optim: general-purpose optimization. Nelder-Mead simplex algorithm.
  cat("Subject",x,"estimated \n")
  return(out)
}

library(parallelsugar) # parallel for Windows # source: https://github.com/nathanvan/parallelsugar
library(fOptions)

# P_BayesPMU <- mclapply(as.list(unique(alldat$id2)),P_kalman_optfun,llfun=P_BayesPMU_mLL,lower=c(.001,.001,.001,.001),upper=c(500,500,2,5),mc.preschedule=FALSE,mc.cores=4)
# P_BayesPI <- mclapply(as.list(unique(alldat$id2)),P_kalman_optfun,llfun=P_BayesPI_mLL,lower=c(.001,.001,.001,.001),upper=c(500,500,2,5),mc.preschedule=FALSE,mc.cores=4)
# L_BayesPI <- mclapply(as.list(unique(alldat$id2)),L_kalman_optfun,llfun=L_BayesPI_mLL,lower=c(.001,.001),upper=c(500,500),mc.preschedule=FALSE,mc.cores=4)

# P_BayesUCB <- mclapply(as.list(unique(alldat$id2)),P_kalman_optfun,llfun=P_BayesUCB_mLL,lower=c(.001,.001,.001,.001),upper=c(500,500,2,5),mc.preschedule=FALSE,mc.cores=4)
L_BayesUCB <- mclapply(as.list(unique(alldat$id2)),L_kalman_optfun,llfun=L_BayesUCB_mLL,lower=c(.001,.001,.001),upper=c(10,10,10),mc.preschedule=FALSE,mc.cores=4)

# P_BayesEI <- mclapply(as.list(unique(alldat$id2)),P_kalman_optfun,llfun=P_BayesEI_mLL,lower=c(.001,.001,.001,.001),upper=c(500,500,2,5),mc.preschedule=FALSE,mc.cores=4)
# L_BayesEI <- mclapply(as.list(unique(alldat$id2)),L_kalman_optfun,llfun=L_BayesEI_mLL,lower=c(.001,.001),upper=c(500,500),mc.preschedule=FALSE,mc.cores=4)


time_finished <- sprintf('Ended at: %s',Sys.time())