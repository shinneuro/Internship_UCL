# for models : L_DecaySMf, P_BayesUCB, P_DecaySMf, L_BayesPI, P_BayesPMU, P_DeltaSMf

# for uncertainty(n_BIC) : P_BayesUCB(12), P_BayesPMU(23), P_BayesSMf(27)
real_dat <- read.csv("Speekenbrink_Konstantinidis_2016_Exp1_anonymous.csv")
alldat <- subset(real_dat,id2 != 1) # excluding first participant. reindexing needed.
load("__final_models__.RData")

Ext <- vector("list")
Uct <- vector("list")

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
  return(list(E,S))
}


# L_DecaySMf

# for (i in 1:79) {
#   # par: eta, theta
#   seq <- subset(alldat,id2==(i+1)) # 2 to 80
#   temp <- L_DecaySMf[[i]]$par
#   eta <- 1/(1+exp(-temp[1]))
#   ut <- LinearUtil(seq$payoff)
#   E <- decayLearning(eta,seq$deck,ut)
#   Ext$L_DecaySMf[[i]] <- E
# }







# P_DecaySMf

# for (i in 1:79) {
#   # par: eta, theta, alpha, lambda
#   seq <- subset(alldat,id2==(i+1)) # 2 to 80
#   temp <- P_DecaySMf[[i]]$par
#   eta <- 1/(1+exp(-temp[1]))
#   theta <- exp(temp[2])
#   alpha <- exp(temp[3])
#   lambda <- exp(temp[4])
#   ut <- ProspectUtil(alpha,lambda,seq$payoff)
#   E <- deltaLearning(eta,seq$deck,ut)
#   Ext$P_DecaySMf[[i]] <- E
# }



# L_BayesPI

# for (i in 1:79) {
#   # par: sg_zeta, sg_epsilon
#   seq <- subset(alldat,id2==(i+1)) # 2 to 80
#   temp <- exp(L_BayesPI[[i]]$par)
#   ut <- LinearUtil(seq$payoff)
#   res <- bayesUpdate(temp[1],temp[2],seq$id2,seq$trial,seq$deck,ut)
#   Ext$L_BayesPI[[i]] <- res[[1]]
#   Uct$L_BayesPI[[i]] <- res[[2]]
# }






# L_BayesPMU
# 
# for (i in 1:79) {
#   # par: sg_zeta, sg_epsilon
#   seq <- subset(alldat,id2==(i+1)) # 2 to 80
#   temp <- exp(L_BayesPMU[[i]]$par)
#   ut <- LinearUtil(seq$payoff)
#   res <- bayesUpdate(temp[1],temp[2],seq$id2,seq$trial,seq$deck,ut)
#   Ext$L_BayesPMU[[i]] <- res[[1]]
#   Uct$L_BayesPMU[[i]] <- res[[2]]
# }


# P_DeltaSMf
# 
# for (i in 1:79) {
#   # par: eta, theta, alpha, lambda
#   seq <- subset(alldat,id2==(i+1)) # 2 to 80
#   temp <- P_DeltaSMf[[i]]$par
#   eta <- 1/(1+exp(-temp[1]))
#   theta <- exp(temp[2])
#   alpha <- exp(temp[3])
#   lambda <- exp(temp[4])
#   ut <- ProspectUtil(alpha,lambda,seq$payoff)
#   E <- deltaLearning(eta,seq$deck,ut)
#   Ext$P_DeltaSMf[[i]] <- E
# }

# P_BayesUCB 

for (i in 1:79) {
  # pars: sg_zeta, sg_epsilon, theta, alpha, lambda
  seq <- subset(alldat,id2==(i+1)) # 2 to 80
  temp <- exp(P_BayesUCB[[i]]$par)
  ut <- ProspectUtil(temp[4],temp[5],seq$payoff)
  res <- bayesUpdate(temp[1],temp[2],seq$id2,seq$trial,seq$deck,ut)
  Ext$P_BayesUCB[[i]] <- res[[1]]
  Uct$P_BayesUCB[[i]] <- res[[2]]
}


# P_BayesPMU

for (i in 1:79) {
  # par: sg_zeta, sg_epsilon, alpha, lambda
  seq <- subset(alldat,id2==(i+1)) # 2 to 80
  temp <- exp(P_BayesPMU[[i]]$par)
  ut <- ProspectUtil(temp[3],temp[4],seq$payoff)
  res <- bayesUpdate(temp[1],temp[2],seq$id2,seq$trial,seq$deck,ut)
  Ext$P_BayesPMU[[i]] <- res[[1]]
  Uct$P_BayesPMU[[i]] <- res[[2]]
}


# P_BayesSMf 

for (i in 1:79) {
  # pars: sg_zeta, sg_epsilon, theta, alpha, lambda
  seq <- subset(alldat,id2==(i+1)) # 2 to 80
  temp <- exp(P_BayesSMf[[i]]$par)
  ut <- ProspectUtil(temp[4],temp[5],seq$payoff)
  res <- bayesUpdate(temp[1],temp[2],seq$id2,seq$trial,seq$deck,ut)
  Ext$P_BayesSMf[[i]] <- res[[1]]
  Uct$P_BayesSMf[[i]] <- res[[2]]
}

save(Ext,Uct,file="__Ext_and_Uct_with_uncertainty__.RData")



# Analysis ----------------------------------------------------------------

# Expectancy
# Uncertainty

# compare between maxExt and choice
# Choice with Maximum Ext: Exploitative, In other case: Exploratory
# Exploitative vs. Exploratory




# for people with PMU -----------------------------------------------------

# participants #:
PMU_selected <- vector()
for (i in 1:79) {
  if(selected[i,2] == "P_BayesPMU"){
    PMU_selected <- c(PMU_selected,i)
  }
}

Ext_PMU <- vector("list")
Uct_PMU <- vector("list")
j <- 1
for (i in PMU_selected) {
  Ext_PMU[[j]] <- Ext$P_BayesPMU[[i]]
  Uct_PMU[[j]] <- Uct$P_BayesPMU[[i]]
  j <- j+1
}

RU_PMU <- vector("list")

for (i in 1:length(PMU_selected)) {
  RU_PMU[[i]] <- Uct_PMU[[i]]/rowSums(Uct_PMU[[i]])
}


# Only for the model that People have selected the most.

consist_Ext <- matrix(0,nrow=200,ncol=length(PMU_selected))# to store whether the choice was exploitative or exploratory

for (i in 1:length(PMU_selected)) {
  for (j in 1:200) {
    if(which.max(Ext_PMU[[i]][j,])==alldat$deck[200*(PMU_selected[i]-1)+j]){
      consist_Ext[j,i] <- "max"
    } else{
      consist_Ext[j,i] <- "xplr"
    }
  }
}

choice_ratio_max <- vector(length=length(PMU_selected))


for (i in 1:length(PMU_selected)) {
  max  <- 0
  xplr <- 0
  for (j in 1:200){
    if (consist_Ext[j,i]=="max") {
      max <- max+1
    } else{
      xplr <- xplr+1
    }
  }
  choice_ratio_max[i] <- max/(max+xplr)
}

# RU_PMU
# histogram for uncertainty 

choice_uncertain <- vector(length=4)
for (i in 1:length(PMU_selected)) {
  for (j in 1:200) {
    choice <- alldat$deck[200*(PMU_selected[i]-1)+j]
    
  }
  
  
}

RT_PMU <- matrix(0,nrow=200,ncol=length(PMU_selected))
for (i in 1:length(PMU_selected)) {
  for (j in 1:200) {
    RT_PMU[j,i] <- alldat$rt[200*(PMU_selected[i]-1)+j]
  }
}


save(Ext_PMU,Uct_PMU,RU_PMU,RT_PMU,choice_ratio_max,file="__Ext_and_Uct_PMU__.RData")

# Choice with Maximum RU: Exploratory, Minimum RU: Exploitative?
# How can I relate those things?

# categorize people?
# 1. When people choose the option with low uncertainty
# 2. When people choose highly uncertain choice
# criterion for high and low uncertaintiy?

# follow uncertainty of each option.

# Expectancy vs. Uncertainty.

# How Human choose and How does it related with the choice?

# Correlate RT with RU

# Distribution over RT
# two groups that were separated by choice strategy (Exploitative / Exploratory)
# we have RU over RT



