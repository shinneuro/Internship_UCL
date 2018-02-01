# Initialization ----------------------------------------------------------

real_dat <- read.csv("Speekenbrink_Konstantinidis_2016_Exp1_anonymous.csv")
alldat <- subset(real_dat,id2 != 1) # excluding first participant. reindexing needed.

# Utility of the Reward  ------------------------------------------------------------------

LinearUtil <- function(R) {
  R
}

ProspectUtil <- function(alpha,lambda,R) {
  ut <- vector("double",length=length(R))
  ut[R < 0] <- -lambda*abs(R[R<0])^alpha # only for the elements in vector R that are negative.
  ut[R >= 0] <- R[R >= 0 ]^alpha
  ut
}

# Learning ----------------------------------------------------------------

deltaLearning <- function(eta,choice,util) {
  n <- length(choice)
  E <- matrix(0.0,nrow=n,ncol=4)
  for(i in 1:(n-1)) {
    E[i+1,] <- E[i,] + eta*(choice[i]==1:4)*(util[i] - E[i,])
  }
  return(E)
}

decayLearning <- function(eta, choice, util) {
 n <- length(choice) 
 E <- matrix(0.0, nrow=n, ncol=4)
 for (i in 1:(n-1)) {
    E[i+1,] <- eta*E[i,] + (choice[i]==1:4)*util[i]
 }
 return(E)
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
  return(E)
}


# Choice ------------------------------------------------------------------

softMax <- function(theta,beta=0,E) {
  # theta is consistency
  # beta a (fixed) bias/exploration "bonus"
  p <- exp(theta*E + beta)
  p <- p/rowSums(p) # make sure is normalized
}


# Delta / SMf -------------------------------------------------------------

P_DeltaSMf_mLL <- function(par,data) {
  eta <- 1/(1+exp(-par[1]))
  theta <- exp(par[2])
  alpha <- exp(par[3])
  lambda <- exp(par[4])
  ut <- ProspectUtil(alpha,lambda,data$payoff)
  E <- deltaLearning(eta,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_DeltaSMf_mLL <- function(par,data) {
  eta <- 1/(1+exp(-par[1]))
  theta <- exp(par[2])
  #alpha <- exp(par[3])
  #lambda <- exp(par[4])
  ut <- LinearUtil(data$payoff)
  E <- deltaLearning(eta,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}


# Decay / SMf -------------------------------------------------------------

P_DecaySMf_mLL <- function(par,data){
  eta <- 1/(1+exp(-par[1]))
  theta <- exp(par[2])
  alpha <- exp(par[3])
  lambda <- exp(par[4])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  E <- decayLearning(eta,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_DecaySMf_mLL <- function(par,data){
  eta <- 1/(1+exp(-par[1]))
  theta <- exp(par[2])
  #alpha <- exp(par[3])
  #lambda <- exp(par[4])
  ut <- LinearUtil(data$payoff)
  E <- decayLearning(eta,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}


# Bayesian / SMf ----------------------------------------------------------

P_BayesSMf_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  theta <- exp(par[3])
  alpha <- exp(par[4])
  lambda <- exp(par[5])
  ut <- ProspectUtil(alpha, lambda, data$payoff)
  E <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

L_BayesSMf_mLL <- function(par,data){
  sg2_zeta <- exp(par[1])
  sg2_epsilon <- exp(par[2])
  theta <- exp(par[3])
  #alpha <- exp(par[4])
  #lambda <- exp(par[5])
  ut <- LinearUtil(data$payoff)
  E <- bayesUpdate(sg2_zeta,sg2_epsilon,data$id2,data$trial,data$deck,ut)
  p <- softMax(theta,beta=0,E)
  -2*sum(log(p[cbind(1:nrow(data),data$deck)]))
}

# Optimization ------------------------------------------------------------

optfun <- function(x,llfun,lower,upper) { # for each participant, parameter estimation
  nstart <- 1000
  dat <- subset(alldat,id2==x)
  pars <- t(runif.sobol(nstart,length(lower),1,0,123345)) # t: matrix transpose # runif: generates random deviates # runif.sobol(n, dimension, init, scrambling, seed): Uniform scrambled Sobol sequence
  pars <- t(lower + (upper - lower)*pars)
  pars[,1] <- log(pars[,1]/(1-pars[,1]))
  if(upper[2] >= 1 & lower[2] > 0) pars[,2] <- log(pars[,2]) # theta
  if(upper[2] < 1 & lower[2] > 0) pars[,2] <- log(pars[,2]/(1-log(pars[,2]))) # theta
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

bayes_optfun <- function(x,llfun,lower,upper) { # for each participant, parameter estimation
  nstart <- 1000
  dat <- subset(alldat,id2==x)
  
  # pars: sg_zeta, sg_epsilon, theta, alpha, lambda
  pars <- t(runif.sobol(nstart,length(lower),1,0,123345)) # t: matrix transpose # runif: generates random deviates # runif.sobol(n, dimension, init, scrambling, seed): Uniform scrambled Sobol sequence
  pars <- t(lower + (upper - lower)*pars)
  # pars[,1:2]: sg_zeta, sg_epsilon
  pars[,1:2] <- log(pars[,1:2])
  # pars[,3]: theta
  if(upper[3] >= 1 & lower[3] > 0) pars[,3] <- log(pars[,3])
  if(upper[3] < 1 & lower[3] > 0) pars[,3] <- log(pars[,3]/(1-log(pars[,3])))
  # pars[,4:5]: alpha, lambda
  pars[,4:ncol(pars)] <- log(pars[,4:ncol(pars)])
  
  tll <- rep(NA,nstart)
  for(i in 1:nstart) {
    tll[i] <-  llfun(par=pars[i,],data=dat)
  }
  par <- pars[which.min(tll),]
  out <- optim(par=par,fn=llfun,data=dat,control=list(maxit=1000)) # optim: general-purpose optimization. Nelder-Mead simplex algorithm.
  cat("Subject",x,"estimated \n")
  return(out)
}


# Result ------------------------------------------------------------------

library(parallelsugar) # parallel for Windows # source: https://github.com/nathanvan/parallelsugar
library(fOptions)
P_DeltaSMf <- mclapply(as.list(unique(alldat$id2)),optfun,llfun=P_DeltaSMf_mLL,lower=c(.001,.001,.001,.001),upper=c(.999,5,2,5),mc.preschedule=FALSE,mc.cores=4)
L_DeltaSMf <- mclapply(as.list(unique(alldat$id2)),optfun,llfun=L_DeltaSMf_mLL,lower=c(.001,.001,.001,.001),upper=c(.999,5,2,5),mc.preschedule=FALSE,mc.cores=4)

P_DecaySMf <- mclapply(as.list(unique(alldat$id2)),optfun,llfun=P_DecaySMf_mLL,lower=c(.001,.001,.001,0),upper=c(.999,5,2,5),mc.preschedule=FALSE,mc.cores=4)
L_DecaySMf <- mclapply(as.list(unique(alldat$id2)),optfun,llfun=L_DecaySMf_mLL,lower=c(.001,.001,.001,0),upper=c(.999,5,2,5),mc.preschedule=FALSE,mc.cores=4)

P_BayesSMf <- mclapply(as.list(unique(alldat$id2)),bayes_optfun,llfun=P_BayesSMf_mLL,lower=c(.001,.001,.001,.001,.001),upper=c(500,500,5,2,5),mc.preschedule=FALSE,mc.cores=4)
L_BayesSMf <- mclapply(as.list(unique(alldat$id2)),bayes_optfun,llfun=L_BayesSMf_mLL,lower=c(.001,.001,.001,.001,.001),upper=c(500,500,5,2,5),mc.preschedule=FALSE,mc.cores=4)

save(P_DeltaSMf,L_DeltaSMf, P_DecaySMf, L_DecaySMf,P_BayesSMf,L_BayesSMf,file="MinhoModels_DeltaSMf.RData")